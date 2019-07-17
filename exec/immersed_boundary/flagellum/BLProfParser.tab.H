/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_BLPROFPARSER_TAB_H_INCLUDED
# define YY_YY_BLPROFPARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif
/* "%code requires" blocks.  */
#line 31 "BLProfParserNC.y" /* yacc.c:1909  */


#include <iostream>
#include <vector>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
using std::ostream;

  struct IVec {
    int iv[BL_SPACEDIM];
    int &operator[](int ii) { return iv[ii]; }
    int get(int ii) { return iv[ii]; }
    std::ostream &operator<<(std::ostream &os) {
      os << '(' AMREX_D_TERM(<< iv[0], << ',' << iv[1], << ',' << iv[2]) << ')';
      return os;
    }
  };

void CPIV(amrex::IntVect &bliv, IVec &ivec);


#line 68 "BLProfParser.tab.H" /* yacc.c:1909  */

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NUMINTEGER = 258,
    NUMFLOAT = 259,
    WS = 260,
    ENDL = 261,
    ENDOFFILE = 262,
    DIGIT = 263,
    ALPHA = 264,
    PLUS = 265,
    MINUS = 266,
    SIGN = 267,
    EXPONENT = 268,
    COMMA = 269,
    EQUAL = 270,
    COLON = 271,
    EOS = 272,
    LPAREN = 273,
    RPAREN = 274,
    SLASH = 275,
    LINE = 276,
    LBRACKET = 277,
    RBRACKET = 278,
    DOTS = 279,
    PLUSSIGN = 280,
    MINUSSIGN = 281,
    DQUOTE = 282,
    SQUOTE = 283,
    DELCOMMA = 284,
    DELENDLINE = 285,
    MPITOKEN = 286,
    OMP = 287,
    WORD = 288,
    UNKNOWN = 289,
    TIME = 290,
    DT = 291,
    COMMENT = 292,
    RUNTIME = 293,
    BLPROFVERSION = 294,
    PHFNAME = 295,
    BLPROFPROC = 296,
    BLPROFDATAFILENAME = 297,
    CALCENDTIME = 298,
    COMMPROFVERSION = 299,
    NPROCS = 300,
    COMMSTATSSIZE = 301,
    CPDATAPROC = 302,
    NCOMMSTATS = 303,
    DATAFILE = 304,
    COMMDATAFILENAME = 305,
    SEEKPOS = 306,
    PROCNAME = 307,
    BARRIERNUMBER = 308,
    NAME = 309,
    QNAME = 310,
    INDEX = 311,
    NAMETAG = 312,
    NAMETAGNAMES = 313,
    TIMEMINMAX = 314,
    TIMERTIME = 315,
    REDUCTION = 316,
    TAGRANGE = 317,
    STEP = 318,
    TIMEGL = 319,
    REGRID = 320,
    WLB = 321,
    LEVEL = 322,
    GRIDS = 323,
    CELLS = 324,
    PCTOD = 325,
    FINESTLEVEL = 326,
    MAXLEVEL = 327,
    REFRATIO = 328,
    PROBDOMAIN = 329,
    COMPUTE = 330,
    SERVICE = 331,
    NOUTFILES = 332,
    HEADERFILE = 333,
    COMMHEADERFILENAME = 334,
    CALLSTATSPROFVERSION = 335,
    CSTATSHEADERFILENAME = 336,
    CSTATSDATAFILENAME = 337,
    REGIONNAME = 338,
    CALLSTATSPROC = 339,
    FNAME = 340,
    INCLUDEALL = 341,
    INCLUDENONE = 342,
    NRSS = 343,
    NTRACESTATS = 344,
    POUND = 345,
    POUNDCOMMENT = 346,
    CPU = 347,
    SLOT = 348,
    CAGE = 349,
    CABINET = 350,
    CAB_POSITION = 351,
    CAB_ROW = 352,
    X_COORD = 353,
    Y_COORD = 354,
    Z_COORD = 355,
    PROCESS_SLOTS = 356,
    PROCESS_SLOTS_FREE = 357,
    PROCESSOR_STATUS_UP = 358,
    PROCESSOR_STATUS_DOWN = 359,
    PROCESSOR_TYPE_SERVICE = 360,
    PROCESSOR_TYPE_COMPUTE = 361,
    ALLOC_MODE_BATCH = 362,
    ALLOC_MODE_OTHER = 363,
    PROCESSOR_ID = 364,
    OD_ALLOCATOR_ID = 365,
    NEXT_RED_BLACK_SWITCH = 366,
    PROCESSOR_SPEC = 367,
    SNULL = 368,
    iv3d = 369
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 55 "BLProfParserNC.y" /* yacc.c:1909  */

  long   iValue;
  double fValue;
  char  *cValue;
  IVec   ivValue;
  std::vector<int> *vintptr;

#line 203 "BLProfParser.tab.H" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void *outpptr);

#endif /* !YY_YY_BLPROFPARSER_TAB_H_INCLUDED  */
