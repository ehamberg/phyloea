# Copyright James H. Bullard 2009

# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the Lesser GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

LEXOR		=	newick.l
GRAMMAR		=	newick.y
PARSER		= 	newick

all:	parser

parser: lex.yy.c
	gcc y.tab.c lex.yy.c -o $(PARSER)

lex.yy.c: y.tab.c $(LEXOR)
	flex $(LEXOR)

y.tab.c: $(GRAMMAR)
	bison -y -d $(GRAMMAR)

clean:
	rm y.tab.c lex.yy.c $(PARSER)
