

## More about `grep`

`grep` is a line by line parser of stdin and by default displays matching lines to the regex pattern.
syntax: 
using stdin: cat file | grep pattern
using files: grep pattern file
common options:
c : count the number of occurrences
m # : repeat match # times
R : recursively through directories
o : only print matching part of line
n : print the line number
v : invert match, print non-matching lines
