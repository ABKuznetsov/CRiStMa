"""Immutable, loss-preserving CIF document nodes."""

from __future__ import annotations

from dataclasses import dataclass

from .tokens import CifToken


@dataclass(frozen=True, slots=True)
class CifScalar:
    tag_token: CifToken
    value_token: CifToken

    @property
    def tag(self) -> str:
        return self.tag_token.value

    @property
    def value(self) -> str:
        return self.value_token.value

    @property
    def raw_value(self) -> str:
        return self.value_token.raw


@dataclass(frozen=True, slots=True)
class CifLoop:
    loop_token: CifToken
    tag_tokens: tuple[CifToken, ...]
    row_tokens: tuple[tuple[CifToken, ...], ...]
    incomplete_tokens: tuple[CifToken, ...] = ()

    @property
    def tags(self) -> tuple[str, ...]:
        return tuple(token.value for token in self.tag_tokens)

    @property
    def rows(self) -> tuple[tuple[str, ...], ...]:
        return tuple(tuple(token.value for token in row) for row in self.row_tokens)

    @property
    def incomplete_values(self) -> tuple[str, ...]:
        return tuple(token.value for token in self.incomplete_tokens)

    def column_index(self, tag: str) -> int | None:
        normalized = tag.casefold()
        for index, candidate in enumerate(self.tags):
            if candidate.casefold() == normalized:
                return index
        return None


@dataclass(frozen=True, slots=True)
class CifBlock:
    name: str
    data_token: CifToken
    scalars: tuple[CifScalar, ...] = ()
    loops: tuple[CifLoop, ...] = ()
    comments: tuple[CifToken, ...] = ()
    unparsed_tokens: tuple[CifToken, ...] = ()

    def scalar(self, tag: str) -> CifScalar | None:
        normalized = tag.casefold()
        return next(
            (
                scalar
                for scalar in self.scalars
                if scalar.tag.casefold() == normalized
            ),
            None,
        )

    def loops_with_tag(self, tag: str) -> tuple[CifLoop, ...]:
        normalized = tag.casefold()
        return tuple(
            loop
            for loop in self.loops
            if any(candidate.casefold() == normalized for candidate in loop.tags)
        )


@dataclass(frozen=True, slots=True)
class CifDocument:
    raw_source: str
    blocks: tuple[CifBlock, ...]
    source_name: str | None = None
    comments: tuple[CifToken, ...] = ()
    unparsed_tokens: tuple[CifToken, ...] = ()
