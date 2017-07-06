

## Reintroducing `grep` (GNU regex parser)

As we have seen in session I, `grep` is a line by line parser by default displays matching lines to the pattern of interest that allows the use of regular expressions (regex) in the specified pattern.

**`grep` usage: **

`cat file | grep pattern`

OR

`grep pattern file`

**`grep` common options:**

`c` : count the number of occurrences

`v` : invert match, print non-matching lines

`R` : recursively through directories

`o` : only print matching part of line

`n` : print the line number

***

## Regular expressions (regex)

Used Pattern matching for a certain amount of text

**Single character:** `B`

**Character sets:**

`[a-z]` : any one from a through z 

`[aei]` : either a, e, i

`[0-9]` : any one from 1 through 9

**Non printable characters:**

`\t` : tab

`\r` : carriage return

`\n` : new line (Unix)

`\r\n` : new line (Windows)

`\s` : space

**Special Characters:**

`.` : *match any character (except new line)*

`\` : *make next character literal*

`^` : *matches at the start of the line*

`$` : *matches at the end of line*

`*` : *repeat match*

`?` : *preceding character is optional*

`( )` : *create a capturing group*

`[ ]` : *sequence of characters*

`{ }` : *place bounds, e.g `{1,6}`*

***

## Examples of regex with `grep`

```bash
grep -c bicycle bicycle.txt
grep "bicycle bicycle" bicycle.txt 
grep ^bicycle bicycle.txt
grep ^Bicycle bicycle.txt 
grep yeah$ bicycle.txt
grep [SJ] bicycle.txt
grep ^[SJ] bicycle.txt 
```
***

## `sed`

`sed` takes a stream of stdin and pattern matches and returns the replaced text to stdout ("Think amped-up Windows Find & Replace").

**`sed` usage:** 

`cat file | sed ‘command’`

OR

`sed ‘command’  file`

**`sed` common options:**

`4d` : *delete line 4*

`2,4d` : *delete lines 2-4*

`/here/d` : *delete line matching here*

`/here/,/there/d` : *delete lines matching here to there*

`s/pattern/text/` : *switch text matching pattern*

`s/pattern/text/g` : *switch text matching pattern globally*

`/pattern/a\text` : *append line with text after matching pattern*

`/pattern/c\text` : *change line with text for matching pattern*

```bash
sed '1,2d' bicycle_copy.txt
sed 's/Superman/Batman/' bicycle_copy.txt 
sed 's/bicycle/car/g' bicycle_copy.txt 
sed 's/.icycle/car/g' bicycle_copy.txt
sed 's/bi*/car/g' bicycle_copy.txt
sed 's/bicycle/tri*cycle/g' bicycle_copy.txt | sed 's/tri*cycle/tricycle
sed 's/bicycle/tri*cycle/g' bicycle_copy.txt | sed 's/tri\*cycle/tricycle
sed '/You/a\best' bicycle_copy.txt
```
