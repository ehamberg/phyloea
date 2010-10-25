/**
   Copyright James H. Bullard 2009

   This program is free software: you can redistribute it and/or modify
   it under the terms of the Lesser GNU General Public License as
   published by the Free Software Foundation, either version 3 of the
   License, or (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   
   You should have received a copy of the Lesser GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

%{

#include <stdio.h>
#include <string.h>
#include "newick.h"
#include "y.tab.h"

%}

%%
([a-zA-Z]+|\"[0-9]+\"|\"[a-zA-Z]+\")              {
				       		       yylval.nameVal = strdup(yytext);
						       return WORD; 
						   }

([0-9]+\.[0-9]*|[0-9]*\.[0-9]+|[0-9]+)            {
						       yylval.numberVal = atof(yytext);
						       return NUMBER;
						  }

[-()<>=+*/;{}:,.]				       return *yytext;

[ \t\n]+					       ;       /* ignore whitespace */

.                                     		       yyerror("Unknown character");
%%



