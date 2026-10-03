from .. import QueryAliases, Structure

def complete(
    source: str,
    cursor: int,
    aliases: QueryAliases | None = ...,
    structure: Structure | None = ...,
) -> tuple[int, int, list[tuple[str, str]]]: ...

__all__ = ["complete"]
