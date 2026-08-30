"""Parser from native CIF tokens to a loss-preserving document."""

from __future__ import annotations

from dataclasses import dataclass, field

from cristma.io.diagnostics import Diagnostic, Severity, SourceSpan
from cristma.io.result import ReadResult, SourceInfo

from .document import CifBlock, CifDocument, CifLoop, CifScalar
from .lexer import lex_cif
from .tokens import CifToken, CifTokenKind


@dataclass(slots=True)
class _BlockBuilder:
    name: str
    data_token: CifToken
    scalars: list[CifScalar] = field(default_factory=list)
    loops: list[CifLoop] = field(default_factory=list)
    comments: list[CifToken] = field(default_factory=list)
    unparsed_tokens: list[CifToken] = field(default_factory=list)

    def freeze(self) -> CifBlock:
        return CifBlock(
            name=self.name,
            data_token=self.data_token,
            scalars=tuple(self.scalars),
            loops=tuple(self.loops),
            comments=tuple(self.comments),
            unparsed_tokens=tuple(self.unparsed_tokens),
        )


def _span(tokens: tuple[CifToken, ...]) -> SourceSpan | None:
    if not tokens:
        return None
    return SourceSpan(tokens[0].span.start, tokens[-1].span.end)


def parse_cif(
    source: str,
    source_name: str | None = None,
) -> ReadResult:
    """Parse CIF source without mapping any crystallographic meaning."""

    tokens, lexer_diagnostics = lex_cif(source)
    diagnostics = list(lexer_diagnostics)
    builders: list[_BlockBuilder] = []
    document_comments: list[CifToken] = []
    document_unparsed: list[CifToken] = []
    current: _BlockBuilder | None = None
    index = 0

    while index < len(tokens):
        token = tokens[index]

        if token.kind is CifTokenKind.COMMENT:
            (current.comments if current is not None else document_comments).append(token)
            index += 1
            continue

        if token.kind is CifTokenKind.DATA:
            current = _BlockBuilder(token.value, token)
            builders.append(current)
            index += 1
            continue

        if current is None:
            document_unparsed.append(token)
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.parse.item_outside_block",
                    "CIF data item appears before the first data_ block",
                    token.span,
                )
            )
            if (
                token.kind is CifTokenKind.TAG
                and index + 1 < len(tokens)
                and tokens[index + 1].kind is CifTokenKind.VALUE
            ):
                document_unparsed.append(tokens[index + 1])
                index += 2
            else:
                index += 1
            continue

        if token.kind is CifTokenKind.TAG:
            next_index = index + 1
            while (
                next_index < len(tokens)
                and tokens[next_index].kind is CifTokenKind.COMMENT
            ):
                current.comments.append(tokens[next_index])
                next_index += 1
            if (
                next_index < len(tokens)
                and tokens[next_index].kind is CifTokenKind.VALUE
            ):
                current.scalars.append(CifScalar(token, tokens[next_index]))
                index = next_index + 1
                continue
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.parse.scalar_value_missing",
                    f"CIF scalar {token.value} has no value",
                    token.span,
                )
            )
            current.unparsed_tokens.append(token)
            index = next_index
            continue

        if token.kind is CifTokenKind.LOOP:
            loop_token = token
            index += 1
            tag_tokens: list[CifToken] = []
            while index < len(tokens):
                candidate = tokens[index]
                if candidate.kind is CifTokenKind.COMMENT:
                    current.comments.append(candidate)
                    index += 1
                    continue
                if candidate.kind is CifTokenKind.TAG:
                    tag_tokens.append(candidate)
                    index += 1
                    continue
                break

            if not tag_tokens:
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.parse.loop_tags_missing",
                        "loop_ is not followed by any data names",
                        loop_token.span,
                    )
                )
                current.unparsed_tokens.append(loop_token)
                continue

            value_tokens: list[CifToken] = []
            while index < len(tokens):
                candidate = tokens[index]
                if candidate.kind is CifTokenKind.COMMENT:
                    current.comments.append(candidate)
                    index += 1
                    continue
                if candidate.kind is CifTokenKind.VALUE:
                    value_tokens.append(candidate)
                    index += 1
                    continue
                if candidate.kind is CifTokenKind.STOP:
                    index += 1
                break

            width = len(tag_tokens)
            complete_count = len(value_tokens) // width * width
            row_tokens = tuple(
                tuple(value_tokens[start : start + width])
                for start in range(0, complete_count, width)
            )
            incomplete = tuple(value_tokens[complete_count:])
            if incomplete:
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.parse.loop_width",
                        f"CIF loop row has {len(incomplete)} values for {width} columns",
                        _span(incomplete),
                    )
                )
            current.loops.append(
                CifLoop(
                    loop_token=loop_token,
                    tag_tokens=tuple(tag_tokens),
                    row_tokens=row_tokens,
                    incomplete_tokens=incomplete,
                )
            )
            continue

        current.unparsed_tokens.append(token)
        severity = (
            Severity.WARNING
            if token.kind in {CifTokenKind.GLOBAL, CifTokenKind.SAVE, CifTokenKind.STOP}
            else Severity.ERROR
        )
        diagnostics.append(
            Diagnostic(
                severity,
                "cif.parse.unexpected_token",
                f"Unexpected CIF token {token.raw!r}",
                token.span,
            )
        )
        index += 1

    document = CifDocument(
        raw_source=source,
        source_name=source_name,
        blocks=tuple(builder.freeze() for builder in builders),
        comments=tuple(document_comments),
        unparsed_tokens=tuple(document_unparsed),
    )
    return ReadResult(
        document=document,
        diagnostics=tuple(diagnostics),
        source_info=SourceInfo(name=source_name, format="cif"),
    )
