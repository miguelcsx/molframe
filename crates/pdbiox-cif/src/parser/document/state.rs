/// What the parser is in the middle of.
struct ParseState<'input, S> {
    sink: S,
    findings: Diagnostics,
    /// Item names of the loop being read, in column order.
    loop_tags: Vec<(Box<str>, Box<str>)>,
    /// Where the loop header began.
    loop_span: Option<ByteSpan>,
    /// How many values of the current loop row have been read.
    loop_cursor: usize,
    /// True between `loop_` and its first value.
    collecting_tags: bool,
    loop_started: bool,
    /// The item a bare `tag value` pair is waiting on.
    pending_tag: Option<(Box<str>, Box<str>, ByteSpan)>,
    saw_block: bool,
    input: PhantomData<&'input str>,
}

impl<'input, S> ParseState<'input, S>
where
    S: ValueSink<'input>,
{
    fn new(sink: S) -> Self {
        Self {
            sink,
            findings: Diagnostics::new(),
            loop_tags: Vec::new(),
            loop_span: None,
            loop_cursor: 0,
            collecting_tags: false,
            loop_started: false,
            pending_tag: None,
            saw_block: false,
            input: PhantomData,
        }
    }

    fn token(&mut self, spanned: Spanned<'input>) {
        match spanned.token {
            Token::Block(name) => self.block(name),
            Token::FrameStart(_) | Token::FrameEnd => self.end_loop(),
            Token::Loop => self.begin_loop(spanned.span),
            Token::Tag(tag) => self.tag(tag, spanned.span),
            Token::Value(text, quoting) => self.value(text, quoting, spanned.span),
        }
    }

    fn block(&mut self, name: &str) {
        self.end_loop();
        self.saw_block = true;
        self.sink.block(name);
    }

    fn begin_loop(&mut self, span: ByteSpan) {
        self.end_loop();
        self.collecting_tags = true;
        self.loop_span = Some(span);
    }

    fn tag(&mut self, tag: &str, span: ByteSpan) {
        let (category, item) = split_tag(tag);
        if self.collecting_tags {
            self.loop_tags.push((category.into(), item.into()));
            return;
        }
        self.end_loop();
        self.pending_tag = Some((category.into(), item.into(), span));
    }

    fn value(&mut self, text: &'input str, quoting: Quoting, span: ByteSpan) {
        if !self.saw_block {
            self.findings.push(Diagnostic::new(Code::E1106).at(span));
            self.saw_block = true;
            self.sink.block("");
        }
        if self.collecting_tags {
            self.sink.begin_loop(&self.loop_tags);
            self.loop_started = true;
        }
        self.collecting_tags = false;

        if let Some((category, item, at)) = self.pending_tag.take() {
            self.sink
                .value(ColumnId(0), &category, &item, text, quoting, at);
            self.sink.end_row();
            return;
        }
        if self.loop_tags.is_empty() {
            self.findings.push(Diagnostic::new(Code::E1105).at(span));
            return;
        }
        // The cursor wraps by comparison rather than by remainder. Two integer
        // divisions per value is a poor price for a counter that only ever
        // advances by one, and a coordinate loop has billions of values.
        let Some((category, item)) = self.loop_tags.get(self.loop_cursor) else {
            return;
        };
        let Some(column) = ColumnId::from_position(self.loop_cursor) else {
            self.findings.push(Diagnostic::new(Code::E1903).at(span));
            return;
        };
        self.sink.value(column, category, item, text, quoting, span);
        self.loop_cursor += 1;
        if self.loop_cursor == self.loop_tags.len() {
            self.loop_cursor = 0;
            self.sink.end_row();
        }
    }

    /// Closes the loop being read, checking that its rows are whole.
    fn end_loop(&mut self) {
        if !self.loop_tags.is_empty() && self.loop_cursor != 0 {
            let columns = self.loop_tags.len();
            let short = self.loop_cursor;
            let mut finding = Diagnostic::new(Code::E1103)
                .with_context("items", columns.to_string())
                .with_context("values in the final row", short.to_string());
            if let Some(span) = self.loop_span {
                finding = finding.at(span);
            }
            if let Some((category, _)) = self.loop_tags.first() {
                finding = finding.in_category(category.to_string());
            }
            self.findings.push(finding);
        }
        if self.loop_started {
            self.sink.end_loop();
        }
        self.loop_tags.clear();
        self.loop_cursor = 0;
        self.loop_span = None;
        self.collecting_tags = false;
        self.loop_started = false;
    }

    fn finish(mut self) -> Result<(S::Output, Vec<Diagnostic>), Vec<Diagnostic>> {
        self.end_loop();
        if self.pending_tag.is_some() {
            self.findings.push(
                Diagnostic::new(Code::E1105)
                    .with_message("an item name was not followed by a value"),
            );
        }
        if !self.saw_block {
            self.findings.push(Diagnostic::new(Code::E1106));
            return Err(self.findings.finish());
        }
        Ok((self.sink.finish(), self.findings.finish()))
    }
}
