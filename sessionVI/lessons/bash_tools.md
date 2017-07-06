

## More about `grep` (GNU regex parser)

`grep` is a line by line parser by default displays matching lines to the pattern of interest that allows the use of regular expressions (regex) in the specified pattern.

**Syntax: **

using stdin: `cat file | grep pattern`

using files: `grep pattern file`

**Common options:**

`c` : count the number of occurrences
`v` : invert match, print non-matching lines
`R` : recursively through directories
`o` : only print matching part of line
`n` : print the line number

## Regular expressions (regex)

Used Pattern matching for a certain amount of text

**Single character:** `B`

**Character sets:**

`[a-z]` # any one from a through z 

`[aei]` # either a, e, i

`[0-9]` : any one from 1 through 9

**Non printable characters:**

`\t` : tab
`\r` : carriage return
`\n` : new line (Unix)
`\r\n` : new line (Windows)
`\s` : space

**Special Characters:**

`.` period or dot: match any character (except new line) 

`\` backslash: make next character literal

`^` caret: matches at the start of the line

`$` dollar sign: matches at the end of line

`*` asterisk or star: *repeat match*

`?` question mark: *preceding character is optional*

`( )` parentheses: create a capturing group

`[ ]` square bracket: sequence of characters

`{ }` curly brace: place bounds, e.g `{1,6}`

