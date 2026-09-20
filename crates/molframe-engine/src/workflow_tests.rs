use crate::{Cost, OperationMetadata, Workflow, WorkflowError, WorkflowInputs};
use molframe_core::{ExecutionContext, MemoryBudget, ScratchPolicy};

fn arithmetic_workflow() -> (crate::CompiledWorkflow, crate::Output<i64>) {
    let mut workflow = Workflow::new();
    let left = match workflow.input::<i64>("left", Cost::Borrow) {
        Ok(value) => value,
        Err(error) => panic!("input failed: {error}"),
    };
    let right = match workflow.input::<i64>("right", Cost::Borrow) {
        Ok(value) => value,
        Err(error) => panic!("input failed: {error}"),
    };
    let sum = match workflow.map2(
        left,
        right,
        OperationMetadata::new("add", Cost::Materialize).with_cse_key("add"),
        |left, right, _| Ok(left + right),
    ) {
        Ok(value) => value,
        Err(error) => panic!("operation failed: {error}"),
    };
    let output = match workflow.output("sum", sum) {
        Ok(value) => value,
        Err(error) => panic!("output failed: {error}"),
    };
    let compiled = match workflow.compile() {
        Ok(value) => value,
        Err(error) => panic!("compile failed: {error}"),
    };
    (compiled, output)
}

#[test]
fn compiled_workflow_is_reusable_and_typed() {
    let (compiled, output) = arithmetic_workflow();
    for value in [2, 9] {
        let mut inputs = WorkflowInputs::new();
        inputs.insert("left", value);
        inputs.insert("right", 3_i64);
        let results = match compiled.run(&inputs, &ExecutionContext::default()) {
            Ok(value) => value,
            Err(error) => panic!("run failed: {error}"),
        };
        match results.get(&output) {
            Ok(result) => assert_eq!(*result, value + 3),
            Err(error) => panic!("typed result failed: {error}"),
        }
    }
}

#[test]
fn missing_and_mistyped_inputs_are_rejected() {
    let (compiled, _) = arithmetic_workflow();
    let mut inputs = WorkflowInputs::new();
    inputs.insert("left", 1_u64);
    assert!(matches!(
        compiled.run(&inputs, &ExecutionContext::default()),
        Err(WorkflowError::InputType { .. })
    ));
}

#[test]
fn missing_input_is_reported_by_name() {
    let (compiled, _) = arithmetic_workflow();
    let mut inputs = WorkflowInputs::new();
    inputs.insert("left", 1_i64);
    assert!(matches!(
        compiled.run(&inputs, &ExecutionContext::default()),
        Err(WorkflowError::MissingInput(name)) if name.as_ref() == "right"
    ));
}

#[test]
fn multiple_outputs_keep_their_types_and_names() {
    let mut workflow = Workflow::new();
    let input = match workflow.input::<u64>("value", Cost::Borrow) {
        Ok(value) => value,
        Err(error) => panic!("input failed: {error}"),
    };
    let doubled = match workflow.map(
        input,
        OperationMetadata::new("double", Cost::Materialize),
        |value, _| Ok(value * 2),
    ) {
        Ok(value) => value,
        Err(error) => panic!("operation failed: {error}"),
    };
    let source = match workflow.output("source", input) {
        Ok(value) => value,
        Err(error) => panic!("output failed: {error}"),
    };
    let result = match workflow.output("doubled", doubled) {
        Ok(value) => value,
        Err(error) => panic!("output failed: {error}"),
    };
    let compiled = match workflow.compile() {
        Ok(value) => value,
        Err(error) => panic!("compile failed: {error}"),
    };
    let mut inputs = WorkflowInputs::new();
    inputs.insert("value", 7_u64);
    let outputs = match compiled.run(&inputs, &ExecutionContext::default()) {
        Ok(value) => value,
        Err(error) => panic!("run failed: {error}"),
    };
    match outputs.get(&source) {
        Ok(value) => assert_eq!(*value, 7),
        Err(error) => panic!("source output failed: {error}"),
    }
    match outputs.get(&result) {
        Ok(value) => assert_eq!(*value, 14),
        Err(error) => panic!("result output failed: {error}"),
    }
    assert!(outputs.get_named::<i64>("doubled").is_err());
}

#[test]
fn compile_rejects_cycles_even_in_a_corrupt_graph() {
    let mut workflow = Workflow::new();
    let input = match workflow.input::<u64>("value", Cost::Borrow) {
        Ok(value) => value,
        Err(error) => panic!("input failed: {error}"),
    };
    let node = match workflow.map(
        input,
        OperationMetadata::new("identity", Cost::Borrow),
        |value, _| Ok(*value),
    ) {
        Ok(value) => value,
        Err(error) => panic!("operation failed: {error}"),
    };
    workflow.nodes[input.index].dependencies.push(node.index);
    workflow.nodes[input.index]
        .expected_inputs
        .push(std::any::TypeId::of::<u64>());
    if let Err(error) = workflow.output("value", node) {
        panic!("output failed: {error}");
    }
    assert!(matches!(workflow.compile(), Err(WorkflowError::Cycle)));
}

#[test]
fn compile_eliminates_common_and_dead_nodes() {
    let mut workflow = Workflow::new();
    let input = match workflow.input::<usize>("value", Cost::Borrow) {
        Ok(value) => value,
        Err(error) => panic!("input failed: {error}"),
    };
    let metadata = OperationMetadata::new("double", Cost::Materialize).with_cse_key("double");
    let first = match workflow.map(input, metadata.clone(), |value, _| Ok(value * 2)) {
        Ok(value) => value,
        Err(error) => panic!("operation failed: {error}"),
    };
    let second = match workflow.map(input, metadata, |value, _| Ok(value * 2)) {
        Ok(value) => value,
        Err(error) => panic!("operation failed: {error}"),
    };
    let _dead = match workflow.map(
        input,
        OperationMetadata::new("dead", Cost::Copy),
        |value, _| Ok(value + 1),
    ) {
        Ok(value) => value,
        Err(error) => panic!("operation failed: {error}"),
    };
    if let Err(error) = workflow.output("first", first) {
        panic!("output failed: {error}");
    }
    if let Err(error) = workflow.output("second", second) {
        panic!("output failed: {error}");
    }
    let compiled = match workflow.compile() {
        Ok(value) => value,
        Err(error) => panic!("compile failed: {error}"),
    };
    assert_eq!(compiled.explain().common_subexpressions_eliminated, 1);
    assert_eq!(compiled.explain().dead_nodes_eliminated, 1);
}

#[test]
fn cancellation_and_memory_admission_happen_before_kernels() {
    let mut workflow = Workflow::new();
    let input = match workflow.input::<usize>("value", Cost::Borrow) {
        Ok(value) => value,
        Err(error) => panic!("input failed: {error}"),
    };
    let result = match workflow.map(
        input,
        OperationMetadata::new("large", Cost::Materialize).with_estimated_output_bytes(1024),
        |value, _| Ok(*value),
    ) {
        Ok(value) => value,
        Err(error) => panic!("operation failed: {error}"),
    };
    if let Err(error) = workflow.output("result", result) {
        panic!("output failed: {error}");
    }
    let compiled = match workflow.compile() {
        Ok(value) => value,
        Err(error) => panic!("compile failed: {error}"),
    };
    let budget = match MemoryBudget::new(32) {
        Ok(value) => value,
        Err(error) => panic!("budget failed: {error}"),
    };
    let context = match ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(0))
        .build()
    {
        Ok(value) => value,
        Err(error) => panic!("context failed: {error}"),
    };
    let mut inputs = WorkflowInputs::new();
    inputs.insert("value", 1_usize);
    assert!(matches!(
        compiled.run(&inputs, &context),
        Err(WorkflowError::MemoryRefused { .. })
    ));

    let context = ExecutionContext::default();
    context.cancellation().cancel();
    assert!(matches!(
        compiled.run(&inputs, &context),
        Err(WorkflowError::Cancelled)
    ));
}
