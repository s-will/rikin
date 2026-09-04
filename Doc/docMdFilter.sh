#!/bin/bash
# Filters a markdown file for Doxygen:
#  - converts LaTeX/Pandoc-style math delimiters to Doxygen's \f$ / \f[ \f]
#  - rewrites absolute GitHub URLs pointing at our own docs pages into bare
#    #anchor fragments, so Doxygen resolves them as internal page links
#    instead of navigating away to github.com (kept working outside Doxygen
#    since the *source* files still use the real GitHub URLs; only this
#    filtered copy, fed to Doxygen, is rewritten)
#
# Fenced code blocks (```...```) are passed through completely unmodified,
# so literal '$', '(', ')', '[', ']' in shell/code examples are never
# mistaken for math delimiters.
#
# Written in plain POSIX-compatible awk (no gensub()/GNU-only extensions),
# since the default /usr/bin/awk on Ubuntu (incl. GitHub Actions runners)
# is mawk, not gawk.

awk '
# Replace every non-overlapping occurrence of text delimited by literal
# strings `open` and `close` (search for the *next* close after each open;
# an open with no matching close is left untouched, same as the original
# sed rules it replaces).
function replace_delim(s, open, closer, prefix, suffix,    result, pos, rest, cpos) {
    result = ""
    while (1) {
        pos = index(s, open)
        if (pos == 0) { result = result s; return result }
        rest = substr(s, pos + length(open))
        cpos = index(rest, closer)
        if (cpos == 0) { result = result s; return result }
        result = result substr(s, 1, pos - 1) prefix substr(rest, 1, cpos - 1) suffix
        s = substr(rest, cpos + length(closer))
    }
}

BEGIN { in_code = 0 }
/^```/ { in_code = !in_code; print; next }
in_code { print; next }
{
    line = $0
    # Order matters: the $$...$$ / $...$ conversions must run before the
    # \(...\) / \[ / \] ones, since those introduce new "$" characters
    # into the line that a later dollar-rule would otherwise reprocess.
    line = replace_delim(line, "$$", "$$", "\\f[", "\\f]")
    line = replace_delim(line, "$", "$", "\\f$", "\\f$")
    line = replace_delim(line, "\\(", "\\)", "\\f$", "\\f$")
    gsub(/\\\[/, "\\f[", line)
    gsub(/\\\]/, "\\f]", line)
    gsub("https://github\\.com/s-will/rikin/blob/master/Help\\.md#", "#", line)
    gsub("https://github\\.com/s-will/rikin#readme", "#rnainterkin-rikin", line)
    gsub("https://github\\.com/s-will/rikin#", "#", line)
    print line
}
' "$1"
