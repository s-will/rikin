#!/bin/bash

cat $1 | sed 's/\\(\(.*\)\\)/\\f$\1\\f$/g;s/\\\[/\\f[/g;s/\\\]/\\f]/g;s/$$\(.\+\)$\$/\\f[\1\\f]/g;s/$\([^$]\+\)\$/\\f$\1\\f$/g;s@https://github\.com/s-will/rikin/blob/master/CliReference\.md#@#@g;s@https://github\.com/s-will/rikin#readme@#rnainterkin-rikin@g;s@https://github\.com/s-will/rikin#@#@g'
