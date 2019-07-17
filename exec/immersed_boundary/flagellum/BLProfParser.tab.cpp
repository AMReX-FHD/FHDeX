/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 4 "BLProfParserNC.y" /* yacc.c:339  */

#include <stdio.h>
#include <AMReX_BLProfStats.H>
#include <AMReX_CommProfStats.H>
#include <AMReX_RegionsProfStats.H>


#define theBLPptr ((BLProfStats *) outpptr)

#define CCCOMMENT SLASH(/)
#define SLASH(s) /##s

#define vout0 if(theBLPptr->Verbose() >= 0) cout
#define vout1 if(theBLPptr->Verbose() >= 1) cout
#define vout2 if(theBLPptr->Verbose() >= 2) cout

int yyerror(void *outpptr, const char *s);
extern int yylex();
extern int yylineno;
extern char *yytext;

static int lineCount = 0;


#line 91 "BLProfParser.tab.cpp" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "BLProfParser.tab.H".  */
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
#line 31 "BLProfParserNC.y" /* yacc.c:355  */


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


#line 145 "BLProfParser.tab.cpp" /* yacc.c:355  */

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
#line 55 "BLProfParserNC.y" /* yacc.c:355  */

  long   iValue;
  double fValue;
  char  *cValue;
  IVec   ivValue;
  std::vector<int> *vintptr;

#line 280 "BLProfParser.tab.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void *outpptr);

#endif /* !YY_YY_BLPROFPARSER_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 297 "BLProfParser.tab.cpp" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  100
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   394

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  115
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  11
/* YYNRULES -- Number of rules.  */
#define YYNRULES  82
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  316

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   369

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   106,   106,   110,   117,   120,   126,   129,   132,   135,
     140,   146,   150,   154,   158,   162,   166,   171,   178,   189,
     194,   199,   204,   209,   214,   227,   243,   259,   266,   272,
     277,   283,   288,   293,   297,   302,   307,   312,   317,   322,
     327,   352,   357,   362,   367,   377,   392,   401,   410,   416,
     425,   432,   447,   455,   462,   469,   476,   482,   488,   493,
     535,   539,   545,   549,   555,   559,   566,   569,   598,   599,
     603,   607,   610,   613,   616,   619,   622,   625,   628,   631,
     634,   637,   640
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUMINTEGER", "NUMFLOAT", "WS", "ENDL",
  "ENDOFFILE", "DIGIT", "ALPHA", "PLUS", "MINUS", "SIGN", "EXPONENT",
  "COMMA", "EQUAL", "COLON", "EOS", "LPAREN", "RPAREN", "SLASH", "LINE",
  "LBRACKET", "RBRACKET", "DOTS", "PLUSSIGN", "MINUSSIGN", "DQUOTE",
  "SQUOTE", "DELCOMMA", "DELENDLINE", "MPITOKEN", "OMP", "WORD", "UNKNOWN",
  "TIME", "DT", "COMMENT", "RUNTIME", "BLPROFVERSION", "PHFNAME",
  "BLPROFPROC", "BLPROFDATAFILENAME", "CALCENDTIME", "COMMPROFVERSION",
  "NPROCS", "COMMSTATSSIZE", "CPDATAPROC", "NCOMMSTATS", "DATAFILE",
  "COMMDATAFILENAME", "SEEKPOS", "PROCNAME", "BARRIERNUMBER", "NAME",
  "QNAME", "INDEX", "NAMETAG", "NAMETAGNAMES", "TIMEMINMAX", "TIMERTIME",
  "REDUCTION", "TAGRANGE", "STEP", "TIMEGL", "REGRID", "WLB", "LEVEL",
  "GRIDS", "CELLS", "PCTOD", "FINESTLEVEL", "MAXLEVEL", "REFRATIO",
  "PROBDOMAIN", "COMPUTE", "SERVICE", "NOUTFILES", "HEADERFILE",
  "COMMHEADERFILENAME", "CALLSTATSPROFVERSION", "CSTATSHEADERFILENAME",
  "CSTATSDATAFILENAME", "REGIONNAME", "CALLSTATSPROC", "FNAME",
  "INCLUDEALL", "INCLUDENONE", "NRSS", "NTRACESTATS", "POUND",
  "POUNDCOMMENT", "CPU", "SLOT", "CAGE", "CABINET", "CAB_POSITION",
  "CAB_ROW", "X_COORD", "Y_COORD", "Z_COORD", "PROCESS_SLOTS",
  "PROCESS_SLOTS_FREE", "PROCESSOR_STATUS_UP", "PROCESSOR_STATUS_DOWN",
  "PROCESSOR_TYPE_SERVICE", "PROCESSOR_TYPE_COMPUTE", "ALLOC_MODE_BATCH",
  "ALLOC_MODE_OTHER", "PROCESSOR_ID", "OD_ALLOCATOR_ID",
  "NEXT_RED_BLACK_SWITCH", "PROCESSOR_SPEC", "SNULL", "iv3d", "$accept",
  "infile", "lines", "line", "knownline", "serviceorcompute", "allocmode",
  "processorstatus", "intvect", "words", "word", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369
};
# endif

#define YYPACT_NINF -62

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-62)))

#define YYTABLE_NINF -4

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
       1,   -27,    -9,   -62,   -62,   -49,   -45,   -62,   -62,   -42,
     -62,   -62,   -62,   -62,   -11,    -3,   -62,   -62,   -62,   -62,
     -62,    28,    50,     2,    73,    76,    74,    79,    80,    86,
      88,    91,    40,    92,    93,    95,    98,    89,    90,   105,
     108,   109,   113,   114,   115,   -61,   116,    65,   121,    70,
     -62,   -62,   -62,   111,   127,    99,   -62,   117,   184,   -62,
     -62,   123,   145,   146,   147,   148,   150,   122,   131,   161,
     -62,   -62,   118,   -62,   -62,   -62,   -62,   106,   119,   165,
     -62,   171,   -62,   166,   175,   177,    23,   178,   -62,   -62,
     176,   179,   -62,   -62,   -62,   -62,   186,   107,   193,   198,
     -62,   -62,   -62,   -62,   -62,   -62,   169,   -62,   191,   207,
     213,   -62,   -62,   -62,   224,   225,   -62,   187,   227,   228,
     -62,   -62,   -62,   -62,   168,   217,   218,   167,   233,   -62,
     176,   -62,   234,   -62,   226,   205,   236,   238,   239,   -12,
     211,   194,   195,   -62,   231,   182,   183,   246,   237,   176,
     163,   157,   250,   240,   252,   253,   229,   -62,   230,   254,
     208,    25,   199,   200,   190,   257,   176,   258,   249,   189,
     264,   265,   266,   -62,   -62,   -62,   219,   255,   256,   259,
     260,    46,    -5,   261,   232,   270,   189,   262,   274,   275,
     276,   220,   221,   279,   280,   214,   222,   284,   -62,   -62,
     206,   277,   286,   287,   290,   263,   235,   241,   -62,   -62,
     -62,   -62,   278,   243,   201,   281,   295,   296,   299,   288,
     289,   -62,   302,   291,   292,   -62,   -62,   -62,   305,   306,
     -62,   308,   309,   -62,   -62,   300,   301,   223,   310,   304,
     303,   203,   317,   307,   311,   242,   313,   312,   319,   320,
     314,   315,   321,   244,   318,   322,   323,   330,   316,   325,
     324,   247,   331,   327,   333,   337,   341,   332,   334,   248,
     335,   338,   345,   346,   -62,   340,   267,   342,   349,   344,
     268,   347,   352,   350,   269,   348,   353,   351,   -52,   -62,
     -62,   354,   -50,   -62,   -62,   356,   -38,   -62,   -62,   358,
     216,   359,   357,   361,   251,   362,   363,   364,   271,   365,
     272,   367,   282,   368,   273,   -62
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,    71,    72,     9,     0,     0,    76,    74,    75,
      77,    78,    80,     6,     0,     0,    70,    73,    79,    81,
      82,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      56,    57,    58,     0,     0,     0,     4,     0,     0,    68,
      10,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      16,    17,     0,    19,    20,    21,    23,     0,     0,     0,
      30,     0,    33,     0,     0,     0,     0,     0,    42,    43,
       0,     0,    22,    27,    50,    48,     0,     0,     0,     0,
       1,     5,     7,    71,    75,     8,     0,    69,     0,     0,
       0,    53,    54,    55,     0,     0,    15,     0,     0,     0,
      29,    32,    31,    41,     0,     0,     0,     0,     0,    44,
       0,    49,     0,    52,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    28,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    12,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    11,    14,    18,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    13,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    66,    45,
       0,     0,     0,     0,     0,    24,     0,     0,    39,    38,
      37,    36,     0,     0,     0,     0,     0,     0,    25,     0,
       0,    67,     0,     0,     0,    46,    47,    26,     0,     0,
      51,     0,     0,    35,    34,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    40,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    64,
      65,     0,     0,    60,    61,     0,     0,    62,    63,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    59
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -62,   -62,   -62,   329,   -62,   -62,   -62,   -62,   -59,   204,
     -58
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    54,    55,    56,    57,   295,   299,   291,   129,    58,
      59
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     107,    -2,     1,    60,     2,     3,    64,    61,     4,   197,
      65,     5,     6,    66,   198,     7,     8,     9,    93,   156,
      94,   157,    67,    10,    11,    12,   125,   126,   177,   178,
      68,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    69,    25,    26,    27,    28,    29,   195,
     196,   289,   290,    70,    30,   293,   294,    71,    31,    32,
      33,    34,    35,    36,    37,    38,    62,    63,    39,   297,
     298,   149,    40,    41,    42,    43,    72,    74,    44,    45,
      73,    46,    75,    76,    47,    48,    49,    50,    51,    77,
     166,    78,    52,    53,    79,    80,    81,    82,    83,    -3,
       1,    84,     2,     3,    85,    86,     4,   183,    87,     5,
       6,    88,    89,     7,     8,     9,    90,    91,    92,    95,
      96,    10,    11,    12,    97,    98,    99,   100,   107,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,   108,    25,    26,    27,    28,    29,   102,   109,   110,
     111,   112,    30,   113,   118,   114,    31,    32,    33,    34,
      35,    36,    37,    38,   115,   116,    39,   117,   120,   122,
      40,    41,    42,    43,   119,   121,    44,    45,   123,    46,
     124,   127,    47,    48,    49,    50,    51,   103,     3,   131,
      52,    53,   103,     3,   128,   132,   133,   130,     7,     8,
     104,   134,   135,     7,     8,   104,    10,    11,    12,   136,
     137,    10,    11,    12,   105,   106,   138,    16,    17,    18,
      19,    20,    16,    17,    18,    19,    20,   139,   140,   141,
     142,   143,   144,   145,   146,   147,   148,   150,   152,   153,
     151,   154,   155,   158,   160,   159,   161,   162,   163,   164,
     168,   165,   167,   169,   170,   171,   172,   175,   176,   181,
     182,   184,   173,   174,   185,   179,   180,   187,   188,   189,
     190,   191,   192,   201,   193,   194,   202,   203,   204,   205,
     199,   200,   208,   209,   210,   206,   207,   212,   213,   215,
     216,   214,   211,   217,   222,   223,   218,   221,   225,   226,
     224,   219,   227,   228,   229,   230,   231,   220,   233,   234,
     232,   235,   236,   240,   237,   238,   243,   242,   239,   241,
     244,   245,   250,   251,   254,   301,   258,   249,   252,   253,
     246,   248,   256,   259,   264,   260,   266,   257,   247,   261,
     267,   255,   265,   262,   268,   263,   269,   271,   274,   275,
     270,   272,   279,   273,   276,   283,   287,   278,   280,     0,
     303,   305,   282,   286,   284,   288,   307,   277,   292,   281,
     296,   285,   300,   186,   302,   304,     0,   306,   308,     0,
     310,   312,   309,   314,   101,   311,   315,     0,     0,     0,
       0,     0,     0,     0,   313
};

static const yytype_int16 yycheck[] =
{
      58,     0,     1,    30,     3,     4,    55,    16,     7,    14,
      55,    10,    11,    55,    19,    14,    15,    16,    79,    31,
      81,    33,    33,    22,    23,    24,     3,     4,     3,     4,
      33,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    15,    43,    44,    45,    46,    47,     3,
       4,   103,   104,     3,    53,   105,   106,    55,    57,    58,
      59,    60,    61,    62,    63,    64,    75,    76,    67,   107,
     108,   130,    71,    72,    73,    74,     3,     3,    77,    78,
       4,    80,     3,     3,    83,    84,    85,    86,    87,     3,
     149,     3,    91,    92,     3,    55,     4,     4,     3,     0,
       1,     3,     3,     4,    15,    15,     7,   166,     3,    10,
      11,     3,     3,    14,    15,    16,     3,     3,     3,     3,
      55,    22,    23,    24,     3,    55,    15,     0,   186,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    18,    43,    44,    45,    46,    47,    30,     3,     3,
       3,     3,    53,     3,    48,    33,    57,    58,    59,    60,
      61,    62,    63,    64,    33,     4,    67,    49,     3,     3,
      71,    72,    73,    74,    55,     4,    77,    78,     3,    80,
       3,     3,    83,    84,    85,    86,    87,     3,     4,     3,
      91,    92,     3,     4,    18,    88,     3,    18,    14,    15,
      16,     3,    33,    14,    15,    16,    22,    23,    24,    18,
       3,    22,    23,    24,    30,    31,     3,    33,    34,    35,
      36,    37,    33,    34,    35,    36,    37,     3,     3,    42,
       3,     3,    64,    16,    16,    68,     3,     3,    33,     3,
      14,     3,     3,    32,    49,    51,    15,    65,    65,     3,
      93,    14,    89,     3,    14,     3,     3,     3,    50,    69,
       3,     3,    33,    33,    15,    66,    66,     3,     3,     3,
      51,    16,    16,     3,    15,    15,    14,     3,     3,     3,
      19,    49,     3,     3,    70,    65,    65,     3,    82,     3,
       3,    14,    70,     3,    51,    94,    33,    19,     3,     3,
      19,    66,     3,    15,    15,     3,    15,    66,     3,     3,
      18,     3,     3,     3,    14,    14,   113,    14,    95,    15,
       3,    14,     3,     3,     3,   109,     3,    15,    14,    14,
      19,    18,    14,     3,     3,    19,     3,    15,    96,    14,
       3,    97,    15,    19,     3,    98,    14,    99,     3,     3,
      16,    16,     3,    15,    14,     3,     3,    15,    14,    -1,
       3,   110,    15,    15,    14,    14,     3,   100,    14,   101,
      14,   102,    14,   169,    15,    14,    -1,    15,    14,    -1,
      15,    14,   111,    15,    55,   113,   113,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   112
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     7,    10,    11,    14,    15,    16,
      22,    23,    24,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    43,    44,    45,    46,    47,
      53,    57,    58,    59,    60,    61,    62,    63,    64,    67,
      71,    72,    73,    74,    77,    78,    80,    83,    84,    85,
      86,    87,    91,    92,   116,   117,   118,   119,   124,   125,
      30,    16,    75,    76,    55,    55,    55,    33,    33,    15,
       3,    55,     3,     4,     3,     3,     3,     3,     3,     3,
      55,     4,     4,     3,     3,    15,    15,     3,     3,     3,
       3,     3,     3,    79,    81,     3,    55,     3,    55,    15,
       0,   118,    30,     3,    16,    30,    31,   125,    18,     3,
       3,     3,     3,     3,    33,    33,     4,    49,    48,    55,
       3,     4,     3,     3,     3,     3,     4,     3,    18,   123,
      18,     3,    88,     3,     3,    33,    18,     3,     3,     3,
       3,    42,     3,     3,    64,    16,    16,    68,     3,   123,
       3,    14,    33,     3,     3,     3,    31,    33,    32,    51,
      49,    15,    65,    65,     3,    14,   123,    89,    93,     3,
      14,     3,     3,    33,    33,     3,    50,     3,     4,    66,
      66,    69,     3,   123,     3,    15,   124,     3,     3,     3,
      51,    16,    16,    15,    15,     3,     4,    14,    19,    19,
      49,     3,    14,     3,     3,     3,    65,    65,     3,     3,
      70,    70,     3,    82,    14,     3,     3,     3,    33,    66,
      66,    19,    51,    94,    19,     3,     3,     3,    15,    15,
       3,    15,    18,     3,     3,     3,     3,    14,    14,    95,
       3,    15,    14,   113,     3,    14,    19,    96,    18,    15,
       3,     3,    14,    14,     3,    97,    14,    15,     3,     3,
      19,    14,    19,    98,     3,    15,     3,     3,     3,    14,
      16,    99,    16,    15,     3,     3,    14,   100,    15,     3,
      14,   101,    15,     3,    14,   102,    15,     3,    14,   103,
     104,   122,    14,   105,   106,   120,    14,   107,   108,   121,
      14,   109,    15,     3,    14,   110,    15,     3,    14,   111,
      15,   113,    14,   112,    15,   113
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   115,   116,   116,   117,   117,   118,   118,   118,   118,
     118,   119,   119,   119,   119,   119,   119,   119,   119,   119,
     119,   119,   119,   119,   119,   119,   119,   119,   119,   119,
     119,   119,   119,   119,   119,   119,   119,   119,   119,   119,
     119,   119,   119,   119,   119,   119,   119,   119,   119,   119,
     119,   119,   119,   119,   119,   119,   119,   119,   119,   119,
     120,   120,   121,   121,   122,   122,   123,   123,   124,   124,
     125,   125,   125,   125,   125,   125,   125,   125,   125,   125,
     125,   125,   125
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     1,     1,     2,     1,     2,     2,     1,
       2,     6,     5,     6,     6,     3,     2,     2,     6,     2,
       2,     2,     2,     2,     8,     9,    10,     2,     4,     3,
       2,     3,     3,     2,    11,    11,     8,     8,     8,     8,
      31,     3,     2,     2,     3,     7,    10,    10,     2,     3,
       2,    10,     3,     3,     3,     3,     1,     1,     1,    65,
       1,     1,     1,     1,     1,     1,     5,     7,     1,     2,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (outpptr, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, outpptr); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, void *outpptr)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (outpptr);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, void *outpptr)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, outpptr);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, void *outpptr)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              , outpptr);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, outpptr); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, void *outpptr)
{
  YYUSE (yyvaluep);
  YYUSE (outpptr);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void *outpptr)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 106 "BLProfParserNC.y" /* yacc.c:1646  */
    { /* empty */
          vout0 << "infile0::lineCount = " << lineCount << endl;
          return 0;
        }
#line 1608 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 3:
#line 110 "BLProfParserNC.y" /* yacc.c:1646  */
    {
          vout0 << "infile1::lineCount = " << lineCount << endl;
          return 0;
        }
#line 1617 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 4:
#line 117 "BLProfParserNC.y" /* yacc.c:1646  */
    {
          ++lineCount;
        }
#line 1625 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 5:
#line 120 "BLProfParserNC.y" /* yacc.c:1646  */
    {
          ++lineCount;
        }
#line 1633 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 6:
#line 126 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "aDELENDLINE" << endl;
      }
#line 1641 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 7:
#line 129 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "kDELENDLINE" << endl;
      }
#line 1649 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 8:
#line 132 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wDELENDLINE" << endl;
      }
#line 1657 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 9:
#line 135 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "infile2::lineCount = " << lineCount << endl;
        vout0 << "ENDOFFILE:  lineCount = " << lineCount << endl;
        return 0;
      }
#line 1667 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 10:
#line 140 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "eDELENDLINE:  " << endl;
      }
#line 1675 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 11:
#line 146 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "kMPI:  " << (yyvsp[-2].iValue) << endl;
           }
#line 1683 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 12:
#line 150 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "kMPI:  " << (yyvsp[-1].iValue) << endl;
           }
#line 1691 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 13:
#line 154 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "kMPI:  " << (yyvsp[-1].iValue) << endl;
           }
#line 1699 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 14:
#line 158 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "kOMP:  " << (yyvsp[-2].iValue) << endl;
           }
#line 1707 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 15:
#line 162 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "RUNTIME:  " << (yyvsp[0].fValue) << endl;
           }
#line 1715 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 16:
#line 166 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "BLPROFVERSION = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->SetBLPVersion((yyvsp[0].iValue));
           }
#line 1724 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 17:
#line 172 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "PHFNAME = " << (yyvsp[0].cValue) << endl;
             theBLPptr->AddFunctionName((yyvsp[0].cValue));
             free((yyvsp[0].cValue));
           }
#line 1734 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 18:
#line 181 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "BLPROFPROC = " << (yyvsp[-4].iValue) << endl;
             vout0 << "DATAFILE   = " << (yyvsp[-2].cValue) << endl;
             vout0 << "SEEKPOS    = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->InitBLProfDataBlock((yyvsp[-4].iValue), (yyvsp[-2].cValue), (yyvsp[0].iValue));
	     free((yyvsp[-2].cValue));
           }
#line 1746 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 19:
#line 189 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "CALCENDTIME = " << (yyvsp[0].fValue) << endl;
	     theBLPptr->AddCalcEndTime((yyvsp[0].fValue));
           }
#line 1755 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 20:
#line 194 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "COMMPROFVERSION = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->SetCPVersion((yyvsp[0].iValue));
           }
#line 1764 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 21:
#line 199 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "NPROCS = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->SetNProcs((yyvsp[0].iValue));
           }
#line 1773 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 22:
#line 204 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "NOUTFILES = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->SetNOutFiles((yyvsp[0].iValue));
           }
#line 1782 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 23:
#line 209 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "COMMSTATSSIZE = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->SetCSSize((yyvsp[0].iValue));
           }
#line 1791 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 24:
#line 218 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "CPDATAPROC = " << (yyvsp[-6].iValue) << endl;
             vout0 << "NCOMMSTATS = " << (yyvsp[-4].iValue) << endl;
             vout0 << "DATAFILE   = " << (yyvsp[-2].cValue) << endl;
             vout0 << "SEEKPOS    = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->InitCommDataBlock((yyvsp[-6].iValue), (yyvsp[-4].iValue), (yyvsp[-2].cValue), (yyvsp[0].iValue));
	     free((yyvsp[-2].cValue));
           }
#line 1804 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 25:
#line 232 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "CPDATAPROC = " << (yyvsp[-7].iValue) << endl;
             vout0 << "NCOMMSTATS = " << (yyvsp[-5].iValue) << endl;
             vout0 << "DATAFILE   = " << (yyvsp[-3].cValue) << endl;
             vout0 << "SEEKPOS    = " << (yyvsp[-1].iValue) << endl;
             vout0 << "PROCNAME   = " << (yyvsp[0].cValue) << endl;
	     theBLPptr->InitCommDataBlock((yyvsp[-7].iValue), (yyvsp[-5].iValue), (yyvsp[-3].cValue), (yyvsp[-1].iValue), (yyvsp[0].cValue));
	     free((yyvsp[-3].cValue));
	     free((yyvsp[0].cValue));
           }
#line 1819 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 26:
#line 248 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "CPDATAPROC   = " << (yyvsp[-8].iValue) << endl;
             vout0 << "NCOMMSTATS   = " << (yyvsp[-6].iValue) << endl;
             vout0 << "COMMDATAFILE = " << (yyvsp[-4].cValue) << endl;
             vout0 << "SEEKPOS      = " << (yyvsp[-2].iValue) << endl;
             vout0 << "NODENAME     = " << (yyvsp[-1].cValue) << endl;
             vout0 << "NODENUMBER   = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->InitCommDataBlock((yyvsp[-8].iValue), (yyvsp[-6].iValue), (yyvsp[-4].cValue), (yyvsp[-2].iValue), (yyvsp[-1].cValue), (yyvsp[0].iValue));
	     free((yyvsp[-4].cValue));
           }
#line 1834 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 27:
#line 260 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "HEADERFILE   = " << (yyvsp[0].cValue) << endl;
	     theBLPptr->AddCommHeaderFileName((yyvsp[0].cValue));
	     free((yyvsp[0].cValue));
           }
#line 1844 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 28:
#line 266 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "BARRIERNUMBER = " << (yyvsp[-2].iValue) << "  NAME = " << (yyvsp[-1].cValue) << "  INDEX = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->AddBarrier((yyvsp[-2].iValue), (yyvsp[-1].cValue), (yyvsp[0].iValue));
	     free((yyvsp[-1].cValue));
           }
#line 1854 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 29:
#line 272 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "NAMETAGINDEX = " << (yyvsp[-1].iValue) << "  INDEX = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->AddNameTag((yyvsp[-1].iValue), (yyvsp[0].iValue));
           }
#line 1863 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 30:
#line 277 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "NAMETAGNAMES = " << (yyvsp[0].cValue) << endl;
	     theBLPptr->AddNameTagName((yyvsp[0].cValue));
	     free((yyvsp[0].cValue));
           }
#line 1873 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 31:
#line 283 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "REDUCTION = " << (yyvsp[-1].iValue) << "  INDEX = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->AddReduction((yyvsp[-1].iValue), (yyvsp[0].iValue));
           }
#line 1882 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 32:
#line 288 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "TIMEMINMAX = " << (yyvsp[-1].fValue) << "  " << (yyvsp[0].fValue) << endl;
	     theBLPptr->AddTimeMinMax((yyvsp[-1].fValue), (yyvsp[0].fValue));
           }
#line 1891 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 33:
#line 293 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "TIMERTIME = " << (yyvsp[0].fValue) << endl;
	     theBLPptr->AddTimerTime((yyvsp[0].fValue));
           }
#line 1900 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 34:
#line 297 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "STEP = " << (yyvsp[-8].iValue) << "  TIME = " << (yyvsp[-5].fValue) << "  WLB = " << (yyvsp[0].iValue) << endl;
	     //theBLPptr->AddTimerTime($2);
           }
#line 1909 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 35:
#line 302 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "STEP = " << (yyvsp[-8].iValue) << "  TIME = " << (yyvsp[-5].iValue) << "  WLB = " << (yyvsp[0].iValue) << endl;
	     //theBLPptr->AddTimerTime($2);
           }
#line 1918 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 36:
#line 307 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "L = " << (yyvsp[-6].iValue) << "  G = " << (yyvsp[-5].iValue) << "  C = " << (yyvsp[-3].iValue) << "  %D = " << (yyvsp[-1].fValue) << endl;
	     theBLPptr->AddGridLevel((yyvsp[-6].iValue), (yyvsp[-5].iValue));
           }
#line 1927 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 37:
#line 312 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "L = " << (yyvsp[-6].iValue) << "  G = " << (yyvsp[-5].iValue) << "  C = " << (yyvsp[-3].iValue) << "  %D = " << (yyvsp[-1].iValue) << endl;
	     theBLPptr->AddGridLevel((yyvsp[-6].iValue), (yyvsp[-5].iValue));
           }
#line 1936 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 38:
#line 317 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "TIME = " << (yyvsp[-5].fValue) << "  WLB = " << (yyvsp[0].iValue) << endl;
	     //theBLPptr->AddTimerTime($2);
           }
#line 1945 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 39:
#line 322 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "TIME = " << (yyvsp[-5].iValue) << "  WLB = " << (yyvsp[0].iValue) << endl;
	     //theBLPptr->AddTimerTime($2);
           }
#line 1954 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 40:
#line 336 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "GLLEV = " << (yyvsp[-30].iValue);
	     vout0 << "  Box = ((" << (yyvsp[-26].iValue) << "," << (yyvsp[-24].iValue) << "," << (yyvsp[-22].iValue) << ") (";
	     vout0 << (yyvsp[-19].iValue) << "," << (yyvsp[-17].iValue) << "," << (yyvsp[-15].iValue) << ") (";
	     vout0 << (yyvsp[-12].iValue) << "," << (yyvsp[-10].iValue) << "," << (yyvsp[-8].iValue) << "))";
	     vout0 << "  len = " << (yyvsp[-5].iValue) << " " << (yyvsp[-4].iValue) << " " << (yyvsp[-3].iValue);
	     vout0 << "  proc = " << (yyvsp[0].iValue);
	     vout0 << endl;
	     theBLPptr->AddGrid3D((yyvsp[-30].iValue),                  // level
	                          (yyvsp[-26].iValue), (yyvsp[-24].iValue), (yyvsp[-22].iValue),          // loend
	                          (yyvsp[-19].iValue) , (yyvsp[-17].iValue) , (yyvsp[-15].iValue),     // hiend
				  (yyvsp[-12].iValue) , (yyvsp[-10].iValue) , (yyvsp[-8].iValue),     // centering
				  (yyvsp[-5].iValue) , (yyvsp[-4].iValue) , (yyvsp[-3].iValue),     // npoints
				  (yyvsp[0].iValue));                // proc
           }
#line 1974 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 41:
#line 352 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "TAGRANGE = " << (yyvsp[-1].iValue) << "  " << (yyvsp[0].iValue) << endl;
	     theBLPptr->AddTagRange((yyvsp[-1].iValue), (yyvsp[0].iValue));
           }
#line 1983 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 42:
#line 357 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "FINESTLEVEL = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->AddFinestLevel((yyvsp[0].iValue));
           }
#line 1992 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 43:
#line 362 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "MAXLEVEL = " << (yyvsp[0].iValue) << endl;
	     theBLPptr->AddMaxLevel((yyvsp[0].iValue));
           }
#line 2001 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 44:
#line 367 "BLProfParserNC.y" /* yacc.c:1646  */
    {
#if (BL_SPACEDIM == 2)
#else
             vout0 << "REFRATIO[lev] = " << (yyvsp[-1].iValue) << " >> " << (yyvsp[0].ivValue)[0] << " << " << endl;
	     amrex::IntVect iv;
             CPIV(iv, (yyvsp[0].ivValue));
	     theBLPptr->AddRefRatio((yyvsp[-1].iValue), iv);
#endif
           }
#line 2015 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 45:
#line 377 "BLProfParserNC.y" /* yacc.c:1646  */
    {
#if (BL_SPACEDIM == 2)
	     //cout << "PROBDOMAIN (spacedim = 2)" << endl;
#else
	     amrex::IntVect ivlo, ivhi, ivc;
             CPIV(ivlo, (yyvsp[-3].ivValue));
             CPIV(ivhi, (yyvsp[-2].ivValue));
             CPIV(ivc,  (yyvsp[-1].ivValue));
	     amrex::Box pd(ivlo, ivhi, ivc);
             //vout0 << "PROBDOMAIN[" << $2 << "] = " << $4 << " " << $5
	           //<< "  " << $6 << endl;
	     theBLPptr->AddProbDomain((yyvsp[-5].iValue), pd);
#endif
           }
#line 2034 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 46:
#line 396 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             cout << (yyvsp[-9].iValue) << " COMPUTE " << (yyvsp[-3].iValue) << " " << (yyvsp[-2].iValue) << " " << (yyvsp[-1].iValue) << " " << (yyvsp[0].iValue) << endl;
	     theBLPptr->AddTopoCoord((yyvsp[-9].iValue), (yyvsp[-3].iValue), (yyvsp[-2].iValue), (yyvsp[-1].iValue), (yyvsp[0].iValue));
           }
#line 2043 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 47:
#line 405 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             cout << (yyvsp[-9].iValue) << " SERVICE " << (yyvsp[-3].iValue) << " " << (yyvsp[-2].iValue) << " " << (yyvsp[-1].iValue) << " " << (yyvsp[0].iValue) << endl;
	     theBLPptr->AddTopoCoord((yyvsp[-9].iValue), (yyvsp[-3].iValue), (yyvsp[-2].iValue), (yyvsp[-1].iValue), (yyvsp[0].iValue), true);
           }
#line 2052 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 48:
#line 411 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "CALLSTATSPROFVERSION = " << (yyvsp[0].iValue) << endl;
             theBLPptr->SetCSVersion((yyvsp[0].iValue));
           }
#line 2061 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 49:
#line 417 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "REGIONNAME = " << (yyvsp[-1].cValue) << endl;
             vout1 << "REGIONNUMBER = " << (yyvsp[0].iValue) << endl;
             theBLPptr->AddRegionName((yyvsp[-1].cValue), (yyvsp[0].iValue));
             free((yyvsp[-1].cValue));
           }
#line 2072 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 50:
#line 426 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "CSTATSHEADERFILENAME   = " << (yyvsp[0].cValue) << endl;
             theBLPptr->AddCStatsHeaderFileName((yyvsp[0].cValue));
             free((yyvsp[0].cValue));
           }
#line 2082 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 51:
#line 437 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "CALLSTATSPROC = " << (yyvsp[-8].iValue) << endl;
             vout1 << "NRSS                 = " << (yyvsp[-6].iValue) << endl;
             vout1 << "NTRACESTATS          = " << (yyvsp[-4].iValue) << endl;
             vout1 << "DATAFILE             = " << (yyvsp[-2].cValue) << endl;
             vout1 << "SEEKPOS              = " << (yyvsp[0].iValue) << endl;
             theBLPptr->InitCStatsDataBlock((yyvsp[-8].iValue), (yyvsp[-6].iValue), (yyvsp[-4].iValue), (yyvsp[-2].cValue), (yyvsp[0].iValue));
             free((yyvsp[-2].cValue));
           }
#line 2096 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 52:
#line 448 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "FNAME = " << (yyvsp[-1].cValue) << "  ";
             vout1 << "FNUMBER = " << (yyvsp[0].iValue) << endl;
             theBLPptr->AddFunctionName((yyvsp[-1].cValue), (yyvsp[0].iValue));
             free((yyvsp[-1].cValue));
           }
#line 2107 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 53:
#line 456 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "PLUS = " << (yyvsp[-1].cValue) << "  " << (yyvsp[0].iValue) << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::ON, (yyvsp[-1].cValue), (yyvsp[0].iValue));
             free((yyvsp[-1].cValue));
           }
#line 2117 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 54:
#line 463 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "MINUS = " << (yyvsp[-1].cValue) << "  " << (yyvsp[0].iValue) << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::OFF, (yyvsp[-1].cValue), (yyvsp[0].iValue));
             free((yyvsp[-1].cValue));
           }
#line 2127 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 55:
#line 470 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "COLON = " << (yyvsp[-1].cValue) << "  " << (yyvsp[0].iValue) << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::UNDEFINED, (yyvsp[-1].cValue), (yyvsp[0].iValue));
             free((yyvsp[-1].cValue));
           }
#line 2137 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 56:
#line 477 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "INCLUDEALL" << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::INCLUDEALL);
           }
#line 2146 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 57:
#line 483 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout1 << "INCLUDENONE" << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::INCLUDENONE);
           }
#line 2155 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 58:
#line 489 "BLProfParserNC.y" /* yacc.c:1646  */
    {
             vout0 << "POUNDCOMMENT :: " << endl;
           }
#line 2163 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 59:
#line 511 "BLProfParserNC.y" /* yacc.c:1646  */
    {
	     vout0 << "_here 0:" << endl;
             vout0 << "CPU = " << (yyvsp[-62].iValue) << endl;
             vout0 << "SLOT = " << (yyvsp[-58].iValue) << endl;
             vout0 << "CAGE = " << (yyvsp[-54].iValue) << endl;
             vout0 << "CAB_POSITION = " << (yyvsp[-46].iValue) << endl;
             vout0 << "CAB_ROW = " << (yyvsp[-42].iValue) << endl;
             vout0 << "X_COORD = " << (yyvsp[-38].iValue) << endl;
             vout0 << "Y_COORD = " << (yyvsp[-34].iValue) << endl;
             vout0 << "Z_COORD = " << (yyvsp[-30].iValue) << endl;
             vout0 << "PROCESS_SLOT = " << (yyvsp[-26].iValue) << endl;
             vout0 << "PROCESS_SLOT_FREE = " << (yyvsp[-22].iValue) << endl;
             vout0 << "PROCESSOR_ID = " << (yyvsp[-12].iValue) << endl;
             vout0 << "OD_ALLOCATOR_ID = " << (yyvsp[-8].iValue) << endl;
	     vout0 << endl;
             theBLPptr->AddEdisonPID((yyvsp[-38].iValue), (yyvsp[-34].iValue), (yyvsp[-30].iValue), (yyvsp[-46].iValue), (yyvsp[-42].iValue), (yyvsp[-54].iValue),
	                             (yyvsp[-58].iValue), (yyvsp[-62].iValue), (yyvsp[-12].iValue));
           }
#line 2186 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 60:
#line 536 "BLProfParserNC.y" /* yacc.c:1646  */
    {
                    vout1 << "PROCESSOR_TYPE_SERVICE" << endl;
		  }
#line 2194 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 61:
#line 540 "BLProfParserNC.y" /* yacc.c:1646  */
    {
                    vout1 << "PROCESSOR_TYPE_COMPUTE" << endl;
		  }
#line 2202 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 62:
#line 546 "BLProfParserNC.y" /* yacc.c:1646  */
    {
                    vout1 << "ALLOC_MODE_BATCH" << endl;
		  }
#line 2210 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 63:
#line 550 "BLProfParserNC.y" /* yacc.c:1646  */
    {
                    vout1 << "ALLOC_MODE_OTHER" << endl;
		  }
#line 2218 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 64:
#line 556 "BLProfParserNC.y" /* yacc.c:1646  */
    {
                    vout1 << "PROCESSOR_STATUS_UP" << endl;
		  }
#line 2226 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 65:
#line 560 "BLProfParserNC.y" /* yacc.c:1646  */
    {
                    vout1 << "PROCESSOR_STATUS_DOWN" << endl;
		  }
#line 2234 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 66:
#line 567 "BLProfParserNC.y" /* yacc.c:1646  */
    {
	 }
#line 2241 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 67:
#line 570 "BLProfParserNC.y" /* yacc.c:1646  */
    {
#if (BL_SPACEDIM == 2)
#else
	   (yyval.ivValue)[0] = (yyvsp[-5].iValue);
	   (yyval.ivValue)[1] = (yyvsp[-3].iValue);
	   (yyval.ivValue)[2] = (yyvsp[-1].iValue);
	   amrex::IntVect iv;
	   iv[0] = (yyvsp[-5].iValue);
	   iv[1] = (yyvsp[-3].iValue);
	   iv[2] = (yyvsp[-1].iValue);
	   //vout0 << "++++ $$ = (" << $2 << "," << $4 << "," << $6 << ")" << endl;
#endif
	 }
#line 2259 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 68:
#line 598 "BLProfParserNC.y" /* yacc.c:1646  */
    { }
#line 2265 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 69:
#line 599 "BLProfParserNC.y" /* yacc.c:1646  */
    { }
#line 2271 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 70:
#line 603 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2<< "wWORD = " << (yyvsp[0].cValue) << endl;
	free((yyvsp[0].cValue));
      }
#line 2280 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 71:
#line 607 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wNUMINTEGER = " << (yyvsp[0].iValue) << endl;
      }
#line 2288 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 72:
#line 610 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wNUMFLOAT = " <<  (yyvsp[0].fValue) << endl;
      }
#line 2296 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 73:
#line 613 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wUNKNOWN" << endl;
      }
#line 2304 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 74:
#line 616 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wEQUAL" << endl;
      }
#line 2312 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 75:
#line 619 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wCOLON" << endl;
      }
#line 2320 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 76:
#line 622 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wCOMMA" << endl;
      }
#line 2328 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 77:
#line 625 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wLBRACKET" << endl;
      }
#line 2336 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 78:
#line 628 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wRBRACKET" << endl;
      }
#line 2344 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 79:
#line 631 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wTIME" << endl;
      }
#line 2352 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 80:
#line 634 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wDOTS" << endl;
      }
#line 2360 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 81:
#line 637 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout2 << "wDT" << endl;
      }
#line 2368 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 82:
#line 640 "BLProfParserNC.y" /* yacc.c:1646  */
    {
        vout1 << "COMMENT:   " << endl;
      }
#line 2376 "BLProfParser.tab.cpp" /* yacc.c:1646  */
    break;


#line 2380 "BLProfParser.tab.cpp" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (outpptr, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (outpptr, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, outpptr);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, outpptr);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (outpptr, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, outpptr);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, outpptr);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 646 "BLProfParserNC.y" /* yacc.c:1906  */


int yyerror(void *outpptr, const char *s) {
  cerr << "*** Unrecognized output at line " << yylineno << "  ::  " << yytext << endl;
  return -1;
}



void CPIV(amrex::IntVect &bliv, IVec &ivec) {
#if (BL_SPACEDIM == 2)
#else
  bliv[0] = ivec[0];
  bliv[1] = ivec[1];
  bliv[2] = ivec[2];
#endif
}

