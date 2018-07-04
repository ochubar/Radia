// From IgorErrors.h.

#ifndef IGOR_ERRORS_H
#define IGOR_ERRORS_H

#define BAD_EXPR_TERM			2001
#define ILLEGAL_POINT_VAL		2002
#define UNKNOWN_DEST_OR_BAD_CMD	2003
#define BAD_EQN_OR_CMD			2004
#define BAD_STR_EXPR			2005
#define STR_TOO_LONG			2006
#define BLD_MAC_ERR				2007
#define X_MAC_ERR				2008
#define MACRO_NESTING_OVERFLOW	2009
#define MACRO_XS_PARAMS			2010
#define MACRO_NO_LPAREN			2011
#define MAC_END_TOO_SOON		2012
#define MAC_PARM_NOT_DEFINED	2013
#define MAC_BAD_PARM_NAME		2014
#define F_BAD_KEYWORD			2015	// unused
#define FIF_ENDIF_MISSING		2016
#define FITTER_LOOP_MISSING		2017
#define F_NO_LPAREN				2018
#define FIF_IF_MISSING			2019
#define FITTER_ITTER_MISSING	2020
#define F_BREAK_NOITTER			2021
#define FDO_DO_MISSING			2022
#define REDUNDANT_PARAM_DEF		2023	// parameters names must be unique
#define NOT_PARAM_NAME			2024	// not a parameter name
#define EXPECT_WAVEPARAM_NAME	2025	// expected a wave name as a function parameter
#define EXPECT_KW_OR_OBJ_NAME	2026	// expected a keyword or an object name
#define EXPECTED_EQU			2027	// expected assignment operator: =, +=, -=, *= or /=
#define NOT_ALLOWED_IN_USER_FCTN 2028	// this function or operation is not available in user functions
#define UNUSED2029				2029	// "unused"
#define EXPECT_COLON_SEP		2030	// "expected : separator"
#define EXPECT_SUBTYPE			2031	// "expected subtype"
#define EXPECTED_EQU_ONLY		2032	// "expected '='"
#define FUNC_END_TOO_SOON		2033	// "unexpected end of function definition"
#define COMP_FUNC_ERROR			2034	// "Function compilation error"
#define COMP_MENU_ERROR			2035	// "Menu compilation error"
#define NO_MACROS_IN_FUNCTIONS	2036	// "Sorry, you can't invoke a macro . . ."
#define CMD_LINE_TOO_LONG		2037	// "The line is too long. Igor command lines are limited to 400 characters."
#define NO_PROMPT_IN_FUNCTIONS	2038	// "Sorry, you can't use Prompt in a function..."

#define NOMEM 1							//  out of memory
#define NOWAV 2							// expected wave name
#define SYNERR 3						// syntax error 
#define NOCOMMA 4						//  expected comma
#define BADAXIS 5						// expected name of active axis
#define BADUNITS 6						//  expected axis units
#define BADNUM 7						//  expected number
#define NOGRAF 8						// there are no graphs
#define BADFLG 9						// unknown flag 
#define BADNPNTS 10						// "the number points in a wave must be between 0 and 20 million."
#define NOEQUALS 11						// missing equal sign in flag (/N=val)
#define NOTPOW2 12						// wave length must be power of 2 for this operation
#define CANT_OVERWRITE_FILE 13			// "file of wrong type -- can not be overwritten"
#define NO_NLFUNC 14					// expected fitting function
#define PNTS_INCOMPATIBLE 15			//  incompatible waveform size
#define NT_INCOMPATIBLE 16				//  incompatible number types
#define NT_FNOT_AVAIL 17				//  can't do function on this number type
#define BADRANGE 18						//  insufficient range
#define NOREALIFFT 19					// can't do IFFT on real data 
#define WNAME_SYNERR 20					// expected wave name
#define NAME_TOO_LONG 21				// name or string too long
#define NO_TERM_QUOTE 22				// expected terminating quote
#define BAD_NAME 23						// ill-formed name
#define NO_SVAR_OR_FCTN 24				// expected string variable or string function
#define NAME_VAR_CONFLICT 25			// name already exists as a variable
#define NAME_CMD_CONFLICT 26			// name already exists as a command
#define NAME_WAV_CONFLICT 27			// name already exists as a waveform
#define NAME_FUNC_CONFLICT 28			// name already exists as a function
#define NAME_MAC_CONFLICT 29			// name already exists as a macro
#define NO_LOG_AXIS 30					// can't do log axis on this axis type
#define NO_SYMFLDR 31					// no symbolic folder of that name
#define NO_SF_DELETE 32					// can't delete special folder
#define TOO_MANY_SFVARS 33				// No more symbolic folders are available.
#define NO_SF_CHANGE 34					// Can't change a special symbolic path.
#define NO_MULT_BIN_SAVE 35				// Can't save multiple waves to a single binary file.
#define BAD_BINARY_FILE 36				// bad or damaged binary file
#define BAD_SPECIAL_HEADER 37			// <<unused>> -- can't find keyword in special file
#define BAD_BIN_VERSION 38				// Bad or incompatible binary version
#define NOKILL_WAVE_IN_USE 39			// can't kill a wave that is part of a graph or table
#define TOO_MANY_PARAMS 40				// Too many parameters on command line.
#define PF_NOCONV 41					// failure to converge during fit
#define BAD_ATTACH_CODE 42				// Improper letter pair for attach code
#define ORN_OFFSET_OUT_OT_RANGE 43 		// Offset out of range
#define NO_ORNS 44						// there are no tags or textboxes
#define NO_NAMED_ORN 45					// there are no tags or textboxes of that name
#define NO_SUCH_FONT 46					// Font not available
#define SSBADCOLUMNNAME 47				// Invalid column name
#define SSNOCOLUMNS 48					// expected column name
#define NO_TARG_WIN 49					// No target window.
#define EXPECT_MORE_DATA 50				// Expecting more data
#define EXPECT_POS_NUM 51 				// expecting a positive non-zero number
#define HIST_APPEND_BAD_DEST 52			// Destination wave must be single precision real for histogram append
#define BAD_MIRROR_AXIS 53				// Mirror axis only available for left or bottom axies when right or top are unused
#define WAVE_NOT_ON_GRAPH 54			// Wave is not on graph
#define NO_REMOVE_CONTROL_WAVE 55 		// can't remove a wave that controls an axis*/
#define MISMATCH_LISTASSIGNMENT 56 		// "mismatch points in wave[a,b]={n,...,n} expression"
#define USER_ABORT 57					// user abort
#define SINGULAR_MATRIX 58				// singular matrix or other numeric error.
#define NO_DATA_IN_FILE 59				// no data was found in file.
#define LINE_TOO_LONG_IN_FILE 60		// file contains a line that is too long.
#define TOOMANYCOLS_IN_FILE 61			// file contains too many columns.
#define WAVE_NAME_TOO_LONG 62			// Wave names can't exceed 31 characters
#define BAD_CHAR_IN_NAME 63				// names must start with a letter and contain letters, digits or '_'
#define TOO_MANY_MAKE_WAVES 64			// More than 100 waves specified in one make command
#define BAD_IGOR_TEXT_FILE 65			// Bad Igor text file. (missing keyword IGOR)
#define BEGIN_BUT_NO_WAVES 66			// Bad Igor text file. (BEGIN keyword but no WAVES)
#define LA_FLAG_PLACEMENT 67			// Flag in wrong place
#define LA_MAX_WAVES 68					// Limit of 12 waves per group exceeded
#define LA_BAD_SYMB_WAVES 69			// Bad symbol on WAVES line
#define LA_NO_NAMES 70					// no wave names found
#define LA_NO_END 71					// missing END or unknown symbol
#define LA_NO_IMAG 72					// bad or missing imaginary number
#define LA_TOO_FEW_PNTS 73				// too little data (less than 2 points)
#define EXPECT_AS_OR_END 74				// Expected 'as', <cr> or ';'
#define DISP_SYNERR 75					// Expected wave name, 'vs', 'as', <cr> or ';'
#define APPEND_SYNERR 76				// Expected wave name, 'vs', <cr> or ';'
#define EXPECT_END 77					// Expected ';' or <cr>
#define EXPECT_WAVE_OR_END 78			// Expected wave name, ';' or <cr>
#define EXPECT_RPAREN 79				// Expected ')' 
#define EXPECT_KEYWORD_OR_NAME 80		// Expected key word or name 
#define DUP_NO_OVERWRITE 81				// Can't overwrite self 
#define EXPECT_LPAREN_OR_EQUAL 82		// Expected '(' or '='
#define UNKNOWN_KEYWORD 83				// Unknown keyword
#define BAD_IRANGE 84					// Expected number between ^2 and ^3
#define PROC_TOO_BIG 85					// The procedure file is getting too big
#define TEXT_TOO_BIG 86					// The window can't hold any more text
#define BAD_FONT_NAME 87				// Font not available
#define NO_PARAMETRIC_X 88				// Can't use a parametric X wave here
#define BAD_CHECKSUM 89					// bad checksum
#define SSMAXDOCSERR 90					// max table documents already opened
#define NO_MORE_PARAMS 91				// no more parameters are expected
#define EXPECTED_XY 92					// expected 'x' or 'y'
#define ONLY_IN_MACRO 93				// this command is only used in macros
#define EXPECTED_VARNAME 94				// expected variable name
#define EXPECTED_STRINGVARNAME 95		// expected string variable name
#define NUM_EXPECTED_NORM 96			// expected comma or end of command (<cr>, ';', or |)
#define NUM_EXPECTED_CRP 97				// expected comma or ')'
#define EXPECT_LPAREN 98				// expected '('
#define EXPECTED_TARG_NAME 99			// expected name of a target window
#define NO_STACK 100					// out of stack space (too much macro recursion)
#define NO_COMPLEX_WAVE 101				// can't use a complex wave here
#define BAD_INDEX 102					// illegal index value
#define NO_LEVEL_FOUND 103				// level crossing not found
#define GRAFWIN_TOO_BIG 104				// graph dimensions too big

#define EXPECT_GRAFWIN 105				// expected name of graph window
#define TOOMANY_OBJECTS 106				// too many graphs, tables or other objects
#define EXPECTED_XOP 107				// expected XOP
#define EXPECTED_XOP_PARAM 108			// expected XOP parameter
#define UNKNOWN_XOP 109					// unknown XOP
#define XOP_OBSOLETE 110				// XOP is incompatible with this version of Igor
#define XOP_HAS_NO_CMD 111				// this XOP has no command associated with it
#define XOP_CANT_INIT 112				// XOP can't initialize itself
#define XOP_INCOMPATIBLE 113			// XOP incompatible with system software or hardware
#define EXPECT_LPAREN_BRACKET 114		// expected '(' or '['
#define INVALID_FORMAT 115				// format string is invalid
#define XOP_NOT_FOUND 116				// can't find file containing XOP
#define EXPECTED_OBJECT 117				// expected graph, table, picture, or textbox name
#define OBJECT_NOT_IN_LAYOUT 118		// object is not in layout
#define NAME_USER_FUNC_CONFLICT 119		// name already exists as user function
#define BAD_HOLD_STRING 120				// hold string length doesn't match number of parameters
#define NO_USER_FCTN 121				// no user defined function of that name exists
#define BAD_FIT_FCTN_NTYPE 122			// functions for fitting must return a single or double precision scalar result
#define BAD_FIT_FCTN_FORMAT 123			// fitting function does not have required form
#define USER_FUC_TOOMANY_PNTS 124		// coefficient wave implys too many fit coefficients (>50)
#define BAD_GROUT 125					// grout must be between 2 and 72 points
#define BAD_TILE_RECT 126				// the tiling window is too small
#define EXPECTED_LAYOUT_OBJECT 127		// expected name of object in layout
#define NO_LAYOUTS 128					// there are no page layouts
#define EBAR_WAVES_MISSING 129			// both positive & negative error bar waves are missing
#define BAD_LAYOUT_EXPAND 130			// magnification must be .25, .5, 1, or 2
#define EXPECTED_LAYOUT 131				// expected name of page layout window
#define NO_PRINT_RECORD 132				// can't open print record (check print driver)
#define ONLY_GRAF 133					// this operation is for graphs only
#define TIME_OUT_READ 134				// timeout while reading data
#define TIME_OUT_WRITE 135				// timeout while writing data

#define BAD_REFNUM 136					// there is no open file with this reference number
#define EXPECTED_TIME 137				// expected time in hh:mm:ss format
#define NOT_TEXT_FILE 138				// non-text files can not be overwritten
#define TOO_MANY_PARAMETERS 139			// too many parameters for this command
#define BAD_WIN_TITLE 140				// expected window title
#define NUM_EXPECTED_RBRACK 141			// expected ']'
#define EXPECT_VERT_AXIS_NAME 142		// "expected vertical axis keyword"
#define EXPECT_HORIZ_AXIS_NAME 143		// "expected horizontal axis keyword"
#define EXPECT_SIZE_KEYWORD 144			// "expected graph size keyword"
#define INDEX_OUT_OF_RANGE 145			// index out of range
#define EXPECTED_FONT 146				// expected font name
#define EXPECTED_MACRO_NAME 147			// expected macro name
#define EXPECT_0_OR_LCBRACE 148			// "expected 0 or '{'"
#define EXPECT_LBRACK 149				// "expected '['

#define EXPECT_LCBRACE 150				// "expected '{'
#define NUM_EXPECTED_RCBRACE 151		// expected '}'

#define EXPECTED_MENU_NAME 152			// expected menu name
#define EXPECTED_MENU_KEYWORD 153		// expected "menu item", SubMenu or End
#define TOO_MANY_SUBMENUS 154			// there are too many submenus
#define CANT_REMOVE_WAVE_TWICE 155		// can't remove the same wave twice
#define FUNCTION_CHANGED_BACK 156		// "A function name cannot be changed in this dialog. It has been changed back to Ò^0Ó."
#define NEED_ANNOTATION_NAME 157		// expected /N=name
#define FLAG_MUST_BE_FIRST 158			// this must be first flag in command
#define NUM_EXPECTED_CRB 159			// expected comma or '}'
#define NAME_WIN_CONFLICT 160			// name already exists as a window name

#define NUM_EXPECTED_CRBRACK 161		// "expected comma or ']'"
#define EXPECTED_GRAPH_OR_TABLE 162		// "expected name of a graph or table"
#define IGOR_OBSOLETE 163				// XOP requires a later version of Igor
#define EXPECT_LINK 164					// "expected link field in options string"
#define EXPECT_LINK_ARG 165				// "expected link argument field in options string"
#define NAME_STRFUNC_CONFLICT 166		// "name already exists as a string function"
#define EXPECT_WAVELIST_LINK 167		// "expected wavelist option link: 'WIN:'"
#define NO_RENAME_BIVARS 168			// "can't rename built-in variables"
#define EXPECT_WAVE_OR_VAR 169			// "expected wave or variable"

#define BAD_OPTIONS_STRING 170			// "inappropriate options string"
#define	DUPLICATE_NAMES 171				// "two or more names are identical"
#define CANT_LOAD_BINARY_FROM_CLIP 172 	// "can't load binary from clipboard"
#define DEMO_ERR 173					// "This is demo..."

#define FIRST_IGOR2_ERR 174				// Start of Igor Pro 2.0 error codes.

#define USR_BKG_ERROR 174				// "Background user function returned error"
#define BKG_NO_TASK_RUNNING 175			// "No background task to control"

#define DASHSPEC_TOO_LONG 176			// "Dash specification too long"

#define BAD_BOUND_EQN 177				// "Wrong format for dependency formula"
#define EXPECT_EQN_DEF 178				// "Expected ':='"
#define NO_PEAK_FOUND 179				// "Peak center not found"
#define EQN_TOO_LONG 180				// "Formula is too long"
#define NO_LOCAL_EQN_VARS 181			// "Local variables can't have dependency formula"

#define ONLY_GRAPH_OR_LAYOUT 182		// this operation is for graphs or layouts only

#define BAD_FPOS_SET 183				// "can't set file mark past end-of-file"
#define NOGRAF_OR_PANEL 184				// "there are no graphs or panels"

#define USING_NULL_STRVAR 185			// "attempt to use a null string variable"
#define FEATURE_NOT_AVAIL 186			// "This feature is not yet available because I don't know how to do it!"

#define INCOMPATIBLE_FLAGS 187			// "incompatible flags"
#define EXPECT_TWO_EXPRS 188			// "expected two expressions"

#define NO_LOCAL_IN_BOUND_EQN 189		// "Can't use local variables in dependency formula"
#define WAVE_TYPE_INCONSISTENT 190		// "Inconsistent type for a wave variable"

#define BAD_FLAG_NUM 191				// "Flag usage is '/N=number or /N=(expression)"

#define BAD_IRANGE_ODD 192				// "Expected odd number between ^2 and ^3"
#define BAD_EXPECTED_RANGE 193			// "expected ^1 between ^2 and ^3" 
#define BAD_ODD_EXPECTED_RANGE 194		// "expected odd ^1 between ^2 and ^3"

#define BAD_XWave 195					// "X data does not match Y data (length or number type)"
#define BAD_WWave 196					// "Weight data does not match Y data (length or number type)"
#define BAD_DWave 197					// "Destination wave does not match Y wave (length or number type)"

#define BAD_LEVELS_SPEC 198				// "missing level specification value(s)"
#define NO_KILL_VAR_OR_STR 199			// "Can't kill variables or strings from a user function"
#define SIZE_OUT_OF_RANGE 200			// "Size out of range"
#define POSITION_OUT_OF_RANGE 201		// "Position out of range"
#define FONTSIZE_OUT_OF_RANGE 202		// "Font size out of range"
#define SHORTER_THAN_COEFS 203			// "Length of wave ^0 can not be less than the coefficient wave"
#define CANT_MIX_CMPLX_WITH_REAL 204	// "Can't combine complex wave ^0 with real wave ^1"
#define CANT_MIX_REAL_WITH_CMPLX 205	// "Can't combine real wave ^0 with complex wave ^1"
#define NAME_ALREADY_IN_USE		 206	// "Name ^0 is already in use"
#define EXPECTED_CONTROL_NAME	 207	// "Expected control name"
#define NO_CHANGE				 208	// "(no change)"
#define EXPECT_MACRO_PROC_FUNC	 209	// "Expected Macro, Proc, or Function keyword"
#define EXPECTED_FUNCTION_NAME	 210	// "Expected function name"
#define TOO_FEW_PARAMETERS		 211	// "Too few parameters"
#define EXPECTED_SUBTYPE		 212	// "Expected subtype"
#define EXPECTED_END_OR_COLON	 213	// "Expected ':' or end of line"
#define WARN_WRONG_CTL_SUBTYPE	 214	// "Warning: wrong subtype for this control, should be Ò : ^0 Ó "
#define WARN_NO_CTLSUBTYPE		 215	// "(Note: optional subtype Ò : ^0 Ó is missing)"	
#define WARN_SUBTYPE_FIXED		 216	// "(Note: optional subtype Ò : ^0 Ó has been added to this procedure)"
#define WARN_PROC_NOT_FOUND		 217	// "Warning: can't find Ò^0Ó in any procedure window"
#define WARN_PROC_EDITED		 218	// "(no change to control, procedure already changed)"
#define NOTE_PROC_EDITED		 219	// "(Note: procedure ^0 has previously been changed)"
#define EXPECTED_NUMEXPR		 220	// "Expected numeric expression"
#define NOQUICKTIME		 		 221	// "QuickTime not present"
#define ERR_MOVIE_ALREADY_OPEN	 222	// "A movie file is already open"
#define ERR_MOVIE_NOT_OPEN		 223	// "No movie file is open"
#define BADSOUNDWAVE			 224	// "Bad sample rate or amplidude for audio wave"
#define NO_WIND_STYLE_MACRO		 225	// "No style macro for this type of window"
#define WRONG_CONTROL			 226	// "^0 is not a ^1 control"
#define EXPECTED_NAME			 227	// "Expected name"
#define RFLAG_NEEDS_CFLAG		 228	// "/C flag must preceed /R flag"
#define ORN_RENAMED				 229	// "(no changes, except ornament already renamed to Ò^0Ó.)"
#define CROSS_AXIS_MISSING		 230	// "Crossing axis not found"
#define CROSS_AXIS_NOT_PERP		 231	// "Crossing axis is not perpendicular"
#define EXPECTED_FIFO			 232	// "expected name of FIFO"
#define FIFO_IN_USE				 233	// "FIFO in use by XOP"
#define NOT_WHILE_FIFO_RUNNING	 234	// "operation not allowed while FIFO is running"
#define NO_FIFO_CHANS			 235	// "no FIFO channels have been defined"
#define FIFO_STOPPED			 236	// "FIFO is not running"
#define NO_SUCH_FIFO_CHAN		 237	// "no FIFO channel of that name"
#define FIFO_OVERFLOW			 238	// "FIFO overflow (disk did not keep up)"
#define WRONG_NUM_CHANS			 239	// "FIFO has a different number of channels"
#define NO_SUCH_ChartCHAN		 240	// "no chart channel of that number"
#define PATH_BUT_NO_FILE		 241	// "/P flag requires file name argument"
#define FILE_NOT_FOUND			 242	// "File not found"
#define EXPECTED_COMPLEX_NUM	 243	// "Expected complex number"
#define EXPECTED_FUNCTION_KEY	 244	// "Expected Function keyword"
#define EXTRA_MACRO_TEXT		 245	// "Extra text in macro, proc, or function"

#define BAD_FILE_TYPE_SPEC		 246	// "A file type specifier must have exactly 4 characters"
#define BAD_FILE_CREATOR_SPEC	 247	// "A file creator specifier must have exactly 4 characters"
#define PATH_TOO_LONG			 248	// "The path to file is too long"
#define FILE_OPEN_READ_ONLY		 249	// "The file ^1 is already open read-only"
#define FILE_OPEN_WRONG_NAME	 250	// "The file ^1 is already open but with the window name ^2"
#define FILE_OPEN_WRONG_TYPE	 251	// "The file ^1 is already open but as a ^2 file"
#define MENU_ITEM_HAS_NO_SUBMENU 252	// "This menu item has no submenu"
#define NO_MACRO_OR_FUNCTION	 253	// "There is no procedure named ^0"
#define CANT_APPEND_TO_THIS_MENU 254	// "Can't add items to this menu"

#define EXPECTED_PICT			 255	// "Expected picture name"
#define CANT_DRAW_HERE			 256	// "Can't draw here"
#define LAYER_NOT_AVAIL			 257	// "Layer not available"
#define NO_DRAW_OBJ				 258	// "No draw object"
#define NO_DOLLAR_HERE			 259	// "Can't use $ here (compiling)"
#define NO_OPEQU_LIST			 260	// "Can't use op= with {n,..,n} assignments"

#define EXPECTED_NOTEBOOK		261		// "Expected name of a notebook window"
#define NB_LOCS_BACKWARDS 		262		// "Invalid notebook selection: the end location is before the start location"
#define NB_STARTLOC_INVALID		263		// "Invalid notebook selection: the start location is out of bounds"
#define NB_ENDLOC_INVALID		264		// "Invalid notebook selection: the end location is out of bounds"

#define EXPECTED_KEYWORD_OR_END 265		// "Expected ',keyword' or end of command (<cr>, ';', or |)"
#define EXPECTED_GRAPHIC_OBJ_NAME 266	// "Expected name of graph, table, layout or picture"
#define NB_BAD_GRAPH_DIMENSIONS 267		// "The graph width and height must be between 50 and 8200"
#define NB_BAD_LAYOUT_RECT		268		// "The layout rectangle is unreasonable"
#define EXPECTED_STRING_EXPR	269		// "Expected string expression"
#define NB_NO_RULER_SPECIFIED	270		// "No ruler specified (use ruler=rulerName or newRuler=rulerName)"
#define TOO_MANY_TABS			271		// "No more than 20 tabs are allowed"
#define BAD_TAB					272		// "Illegal tab value"
#define EXPECTED_SELECTION_KW	273		// "Expected notebook selection keyword"
#define EXPECTED_PATH			274		// "Expected name of a symbolic path"
#define NB_UNKNOWN_SPECIAL_CHAR_TYPE 275 // "Unknown special character code"
#define FUNCTION_TOO_LONG		276		// "The function is too long"
#define BAD_GRAPH_PREFS			277		// "Bad Graph Preferences"
#define BAD_TABLE_PREFS			278		// "Bad Table Preferences"
#define BAD_LAYOUT_PREFS		279		// "Bad Layout Preferences"
#define BAD_PANEL_PREFS			280		// "Bad Panel Preferences"
#define BAD_NOTEBOOK_PREFS		281		// "Bad Notebook Preferences"
#define BAD_PROCEDURE_PREFS		282		// "Bad Procedure Preferences"
#define BAD_COMMAND_PREFS		283		// "Bad Command Window Preferences"
#define BAD_WINDOWS_PREFS		284		// "Bad Windows Preferences"
#define BAD_PALETTE_PREFS		285		// "Bad Color Palette Preferences"
#define BAD_DASHED_PREFS		286		// "Bad Dashed Lines Preferences"
#define BAD_MISC_PREFS			287		// "Bad Miscellaneous Preferences"
#define BAD_HEADER_FOOTER_PREFS	288		// "Bad Header/Footer Preferences"

#define NO_DATA_IN_CLIP			289		// "no data was found in the clipboard"

#define PICT_ALREADY_EXISTS		290		// "a picture by that name already exists"
#define EXPECTED_PICT_NAME		291		// "expected the name of a picture"
#define NO_PICT_IN_CLIP			292		// "there is no picture in the clipboard"
#define RESOURCE_NOT_FOUND		293		// "resource not found"

#define TOOMANY_WLASSIGN_PNTS 	294		// "too many points in wave={n,..,n} expression"
#define NO_KILL_SELF			295		// "A window hook procedure tried to kill its own window"
#define FIT_COEFF_DIFFERENT		296		// "The coefficient wave is not the same number type as the data. Use Redimension."
#define KN_NO_EQN				297		// "Kn variables can't have dependency formulas."
#define STD_AXISNAME_WRONG_EDGE	298		// "Can't use standard axis name on different edge."
#define EXPECT_VERT_AXIS		299		// "Expected vertical axis."
#define EXPECT_HORIZ_AXIS		300		// "Expected horizontal axis."
#define NOT_FIFO_FILE			301		// "Not a FIFO file."
#define BAD_FIFO_FILE_VERS		302		// "Bad FIFO file version."
#define CORRUPT_FIFO_FILE		303		// "Corrupt FIFO file."

#define CANT_KILL_PICT			304		// "can't kill a picture that is used in a graph, layout or panel"
#define EXPECT_WAVE_VAR_STR		305		// "expected name of wave, variable, or string"
#define NAME_PICT_CONFLICT		306		// "name already exists as a picture"
#define BAD_PICT_VERSION		307		// "picture version is not valid"

#define CMD_CAN_NOT_BE_COMPILED	308		// "Sorry, this operation is not allowed in a function. It is allowed in a macro."
#define EXPECTED_TABLE			309		// "expected name of a table"
#define DELIM_RADIX_CONFLICT	310		// "comma can not be both a delimiter and the decimal character"
#define IT_EXPECTED_BEGIN		311		// "expected BEGIN keyword after WAVES in Igor Text file
#define IT_UNKNOWN_KEYWORD		312		// "unknown keyword in Igor Text file"
#define LF_BAD_COLUMN_NAME_LINE	313		// "column names must be before data"
#define LF_BAD_LINE_OR_COL		314		// "line and column numbers must be ³ 0"

#define MENU_HELP_MISPLACED 	315		// "the help must appear on the line after the menu item string expression"
#define MENU_KEYWORD_UNKNOWN 	316		// "the only keyword allowed after 'Menu' is 'dynamic'"

#define NO_IGOR_OBJECT_INFO_PICCOMMENT 317	// "The picture does not contain Igor object information"
#define BAD_IGOR_OBJECT_INFO_PICCOMMENT 318	// "The picture contains Igor object information that this version of Igor does not understand"

#define NO_TABLES				319		// "There are no tables",
#define NO_XFUNC				320		// "No external function of that name exists"
#define NO_USERFUNC_OR_XFUNC	321		// "No user or external function of that name exists"

#define NO_COEFS_WAVE			322		// "You need a coefficients wave for fitting to a user-defined function"
#define BAD_POLY_NTERMS			323		// "Expected a number of polynomial terms from 3 to 20"
#define TOO_MANY_WAVES			324		// "Too many waves - 100 is the maximum"

#define NAME_PATH_CONFLICT 		325		// "Name already exists as symbolic path"
#define RENAME_CONFLICT			326		// "You have already renamed Ò^0Ó to Ò^1Ó."

#define NO_CHANGE_WAVE_INUSE	327		// "Can't change a wave in the middle of a wave assignment."
#define NO_WAVE_X_DEST			328		// "Can't use wave(x)= in a function. Use x2point and wave[p]= instead"
#define EXPECT_MARGIN			329		// "Expected margin keyword (left,right,top or bottom)"
#define NULL_WAVE_OP			330		// "Attempt to operate on a NULL or missing wave"
#define NAME_IS_RESERVED_KW		331		// "Name is a reserved keyword"
#define NOCOMPILE_APPEND		332		// "Can't compile Append. Use AppendToGraph or AppendToLayout or AppendToTable"
#define NOCOMPILE_REMOVE		333		// "Can't compile Remove. Use RemoveFromGraph or RemoveFromLayout or RemoveFromTable"
#define AXISENAB_RANGE_ERR		334		// "Axis enable settings out of range: must be between 0 and 1 and start < stop"
#define NEED_COMPILE			335		// "The procedure window(s) need to be compiled. Perhaps auto-compile is off."
#define NOKILL_OBJ_IN_FCTN		336		// "Can't kill a wave or variable that is part a user function or dependency expression."
#define TAG_FUNCTION_ERROR		337		// "A tag access function is only valid while a tag is being drawn."
#define TRACE_SPECIFED_TWICE	338		// "A trace was specifed twice."

#define WIN_TITLE_BAD_LENGTH	339		// "Window titles can be 1 to 40 characters long"
#define UNKNOWN_LUMP_REC_VERSION 340	// "This version of Igor can't handle this packed experiment file"
#define CANT_UNTIL_CHOOSE_FROM_LIST 341	// "Select an item in the list"
#define XOP_RESOURCES_MISSING	342		// "The XOP file Ò^1Ó is missing required resources. . ."
#define XOP_NEEDS_FPU			343		// "This XOP can't be loaded because it requires a math coprocessor."

#define NUM_BINS_MUST_BE_TWO 344			// "Histogram requires at least two bins"
#define SRC_AND_DEST_MUST_BE_DIFFERENT 345	// "Source and destination waves must be different"
#define BIN_WIDTH_CANT_BE_ZERO 346			// "Histogram bin width can't be zero"
#define BIN_PARAMS_NAN_OR_INF 347			// "Histogram bin start or bin width is a NaN or an INF"

#define RULERS_MUST_BE_DIFFERENT 348		// "The rulers must be different"
#define EXPECTED_GRAPH_OR_LAYOUT 349		// "expected name of a graph or layout"
#define SAFE_SAVE_DISK_FULL 350				// "The save could not be done because the disk is full"
#define DIRECT_XFUNC_CANT_DO_CALLBACK 351	// "XFUNC programmer error: a direct XFUNC is not allowed to do a callback to Igor."

#define INCLUDE_BAD_FILE_TYPE 352			// "#included files must be of type TEXT"
#define BAD_INCLUDE_SPEC 353				// "the #include file specification is bad"
#define INCLUDE_FILE_NOT_FOUND 354			// "include file not found"
#define READONLY_PROCEDURE 355				// "This is a read-only procedure. Change the name to create a new procedure."

#define TU_BAD_PARAGRAPH 356				// "The TU paragraph is out of range."
#define TU_BAD_LOC_CHAR_POS 357				// "The TU location character position is out of range."
#define TU_BAD_LOC_ORDER 358				// "The first TU location is AFTER the second TU location."

#define XOP_RECURSION_ATTEMPTED 359			// "The XOP has attempted recursion"
#define INCLUDED_FILE_OUTDATED 360			// "The included procedure file Ò^1Ó needs to be updated"

#define BUTTON_NEEDS_COMPILED 361			// "You need to compile the procedure windows before you use controls."
#define NO_BUTTON_PROC 362					// "Cannot find control's action procedure Ò^0Ó

#define TOO_LONG_FOR_CMD_LINE 363			// "Too long to fit on command line"
#define CLICK_AUTO_DUMMY 364				// "Click ÒSet to Auto ValuesÓ button"
#define NEED_N_DIGITS 365					// "You need at least ^2 digits to properly label each tick increment."
#define START_EXCEEDS_END 366				// "Start value must be less than end value"
#define PATH_IN_USE 367						// "The symbolic path is in use and can't be killed."

#define XOP_REQUIRES_PPC 368				// "This XOP will run on a PowerMac but not on this 680x0 Mac."
#define CODE_FRAGMENT_LOAD_ERROR 369		// "A system error (%d) occurred while loading the code fragment \"%s\"."

#define XFUNC_BAD_NT 370					// "This XFUNC was compiled with an illegal number type (must be double)."

#define DO_WINDOW_FROM_FUNCTION 371			// "DoWindow/R requires recompiling functions and thus can't be called from a function."

#define FIRST_IGOR3_ERR 372					// Start of Igor Pro 3.0 error codes.

#define NO_TEXTNUMERIC_WAVE_OVERWRITE 372	// "Can't overwrite a text wave with a numeric wave or vice versa."
#define NEED2PNTS_FOR_THIS_OP 373			// "This operation requires a wave with at least two points."
#define NO_TEXT_OP 374						// "This operation does not work on text waves."
#define NODATA_IN_DIM 375					// "There is no data allocated in this dimension."
#define DIMENSION_MISMATCH 376				// "Mismatch between actual and specified number of dimensions."
#define BAD_INCREMENT 377					// "Increment value is less than 1."
#define ZERO_DATA_IN_WAVE 378				// "The wave has zero data allocated."
#define INCONSISTENT_DIMENSIONS 379			// "Inconsistent number of items in a dimension."
#define BAD_DIMENSION_NUMBER 380			// "Dimension number out of range (0-3)."
#define EXPECT_DIM_LABEL 381				// "Expected a dimension item label (literal, not a string)."
#define NOSUCH_DIM_LABEL 382				// "Couldn't find the given dimension item label."
#define NO_PARENT_DATAFOLDER	383			// "Tried to access the parent data folder of the root or of a non-existent data folder."
#define NO_CHILD_DATAFOLDER		384			// "No child data folder of that name exists."
#define CANT_APPEND_DIF_X_TO_CAT 385		// "Can't append different x-wave to category axis."
#define CANT_APPEND_TEXT_X_TO_NONCAT 386	// "Can't append text x-wave to non-category axis."
#define TOO_MANY_SORT_SRC_WAVES	387			// "Too many sort key waves specified. Maximum is 10."
#define CANT_PUT_DF_IN_SELF		388			// "Can't move a data folder into itself or into a child folder."
#define CANT_RENAME_ROOT		389			// "Can't rename the root data folder."
#define EXPECT_DATAFOLDER_NAME		390		// "Expected name of a data folder"
#define NO_GLOBALS_IN_FUNCTIONS		391		// "Must use WAVE, NVAR & SVAR to access global waves, variables & strings with #pragma rtGlobals=1."
#define EXPECT_SQUARE_MAT		392			// "Expected a square matrix."

#define NO_ROOT_DATAFOLDER		393			// "A non-existent data folder was referenced while accessing a child data folder."
#define CANT_FIND_FOLDER		394			// "Can't find data folder."
#define CANT_KILL_PARENTFOLDER_OF_CURRENT 395	// "Can't kill a data folder that contains the current data folder."
#define FOLDER_NAME_EXISTS		396			// "A data folder of that name already exists at this level."
#define EXPECT_COMPILER_DIRECTIVE	397		// "Expected a compiler directive."
#define UNKNOWN_COMPILER_DIRECTIVE	398		// "Unknown compiler directive."
#define EXPECT_PRAGMA_KW		399			// "Expected a pragma keyword."
#define UNKNOWN_PRAGMA_KW		400			// "Unknown pragma keyword."
#define GVAR_TYPE_INCONSISTENT	401			// "Inconsistent type for a global variable reference."

#define OBSOLETE_SCRAP 402					// "The clipboard contents are not compatible with this version of Igor."
#define NUMERIC_WAVE_CANT_HAVE_TEXT_FORMAT 403	// "A numeric wave can not have a text format."
#define TEXT_WAVE_CANT_HAVE_NUMERIC_FORMAT 404	// "A text wave can not have a numeric or date/time format."

#define BAD_CHAR_IN_WAVE_NAME 405			// "A wave name can't contain any of the following: ' \" ; or : "*/
#define BAD_COLORINDEX_WAVE 406				// "Expected a matrix wave containing 3 columns with red, green, and blue values."
#define EXPECT_IMAGE_NAME 407				// "Expected the name of an image in the top graph."
#define EXPECT_MATRIX 408					// "Expected a 2D (matrix) wave."
#define EXPECT_VS_OR_END 409				// "Expected 'vs' keyword or end of command."
#define EXPECT_COLORTABLE_NAME 410			// "Expected name of a color table."

// These errors would normally be generated by a buggy XOP but could also be generated by a buggy Igor.
#define UNKNOWN_WAVE_ACCESS_MODE 411		// "An attempt was made to access a wave using an incompatible access mode."
#define NUMERIC_ACCESS_ON_TEXT_WAVE 412		// "An attempt was made to treat a text wave as if it were a numeric wave."
#define TEXT_ACCESS_ON_NUMERIC_WAVE 413		// "An attempt was made to treat a numeric wave as if it were a text wave."
#define MD_WAVE_BAD_INDEX 414				// "An invalid index was used to access a wave."

#define CONTOUR_EXPECTED_XY_WAVES 415		// "Expected \"vs {xwave,ywave}\"."
#define EXPECT_1DZ_WAVE 416					// "Expected a 1D (single column) contour data (z) wave."
#define EXPECTED_CONTOUR_XYZMATRIX 417		// "Expected a matrix wave containing 3 columns with X, Y, and Z values, or \"zwave vs {xwave,ywave}\"."
#define EXPECTED_1D_XYWAVE 418				// "Expected a 1D (single column) xwave or ywave in \"vs {xwave,ywave}\"."
#define CONTOUR_SHORT_XWAVE 419				// "xwave in \"vs {xwave,ywave}\" has too few rows for contour data rows."
#define CONTOUR_SHORT_YWAVE 420				// "ywave in \"vs {xwave,ywave}\" has too few rows for contour data columns."
#define CONTOUR_MISMATCHED_WAVES 421		// "XYZ waves must all have the same length."

// These errors are used by GetSelection. Prior to Igor Pro 2.5, they were equated to NOMEM because HR forgot to create real error messages for them.
#define EXPECTED_WINTYPE 422						// "Expected a window type keyword: graph, table, layout, notebook or panel."
#define EXPECTED_GRAPH_TABLE_LAYOUT_NOTEBOOK 423	// "Expected one of these window type keywords: graph, table, layout, notebook."
#define EXPECTED_TABLE_WIN 424						// "Expected name of table window"
#define EXPECTED_PANEL 425							// "Expected name of panel window"

#define NO_ROW_COL_LABELS_POSITIONS 426		// "Reading of row or column labels and positions is supported only for delimited-text matrix loading."
#define TOO_FEW_COLUMNS_FOR_MATRIX_LOAD 427	// "There are not enough columns in the file to load the matrix using the specified parameters."

#define LA_ONE_MD_WAVE_PER_BLOCK 428		// "WAVES declaration error: Each multi-dimensional wave must be in its own data block."
#define LA_BAD_DIMENSIONS 429				// "WAVES declaration error: Improper specification of wave dimension sizes."
#define FORMAT_STR_TOO_LONG 430				// "Format string too long."
#define FORMAT_STR_NO_PCT 431				// "Format string missing '%' character."
#define ILLEGAL_CONTOUR_FORMAT 432			// "Format string needs one '%f' or '%g', may have at most one '%d' and one '%s'."

#define MUST_BE_TWO_FREE_DIMENSIONS 433		// "There must be two free dimensions: one vertical (-2) and one horizontal (-3)."
#define ONE_OF_COLOR_COLORINDEX_OR_CTAB 434 // "Expected only one of color=, ctab= or cindex= keywords."
#define SUBDIVISIONS_REQUIRES_XYZ 435		// "The interpolate= keyword applies only to XYZ contour data, yet the named contour data is a Z matrix."
#define ZNULL_REQUIRES_XYZ 436				// "The nullValue= keyword applies only to XYZ contour data, yet the named contour data is a Z matrix."
#define ONLY_ONE_LEVELS_ARG_ALLOWED 437		// "Expected only one autoLevels= or manLevels= keyword."
#define EXPECT_1D_LEVEL_WAVE 438			// "Expected a 1D (single column) wave in \"manLevels=waveName\"."
#define EXPECT_CONTOUR_NAME 439				// "Expected the name of a contour matrix or z wave in the top graph."

#define BAD_OBJECT_TYPE 440					// "Illegal object type."
#define BAD_VAR_INDEX 441					// "Numeric variable index out of range."
#define BAD_STR_INDEX 442					// "String variable index out of range."
#define CONTOURXYZ_TOO_FEW_ROWS 443			// "Need at least 4 rows of X, Y, and Z to generate a contour"

#define INCOMPATIBLE_DIMENSIONING	444		// "Incompatible dimensions."
#define EXPECT_MATRIX_OR_VECTOR		445		// "Expected a matrix or vector and got a higher dimensioned object."
#define OSACompileErr				446		// "Got an error when compiling an OSA Script (AppleScript)."
#define OSAExecuteErr				447		// "Got an error when executing an OSA Script (AppleScript)."
#define EXPECT_STR_EXPR_NO_FUNC		448		// "Expected string expression NOT involving local variables or user functions."
#define CANT_FIND_SCRIPTING_SYSTEM	449		// "Can't connect to the Scripting System (AppleScript). Make sure it was installed correctly."
#define BAD_CHAR_IN_DF_NAME		450			// "A data folder name can't contain any of the following: ' \" ; : or any control chars."
#define AUTO_VAR_CONFLICT		451			// "Conflict creating V_ or S_ local variable. Probably user's NVAR, SVAR or wrong type."
#define FFT_ROWS_EVEN			452			// "rows must be even for FFT on real data"
#define BADSOUNDWAVEFREQ		453			// "Invalid sample rate for sound wave (use SetScale)"
#define OLD_SOUNDMGR			454			// "Sound manager version 3.0 or later is not present"
#define NO_TEXT_YTRACE			455			// "Can't use a text wave as y trace. Use ModifyGraph swapXY to get horizontal bars."
#define NVARorSVARfailed		456			// "Failed to resolve a local reference to a global variable or string (NVAR or SVAR)."
#define NOT_ON_WINDOWS			457			// "Feature not available on Windows."

#define CMD_ENDED_WITHOUT_NAME	458			// "Expected name"
#define CANT_CUT_MIX_OF_ROWS_AND_COLUMNS 459 // "Can't cut a mix of rows from one wave and columns from another."
#define ZNULL_MANUAL_OR_AUTO	460			// "Expected nullValueAuto or nullValue=value, but not both"
#define TOO_MANY_LEVELS			461			// "Too many levels (Max ^3)"

#define WAVE_DIMENSION_OR_VIEW_CONFLICT 462	// "All waves must have the same number of dimensions (at least two dimensions) and same viewed dimensions."
#define TRIANGULATION_REQUIRES_XYZ 463		// "The triangulation= keyword applies only to XYZ contour data, yet the named contour data is a Z matrix."

#define NOT_PACKED_EXP_FILE 464				// "This is not a packed Igor experiment file or it has been corrupted."
#define CONFLICT_DIFFERENT_TYPES 465		// "Can't create '^0' because that name is in use for a different type object."
#define INCOMPATIBLE_DATA_BROWSER 466		// "This version of the Data Browser is not compatible with this version of Igor."
#define SUB_DATA_FOLDER_NOT_FOUND 467		// "The specified sub data folder was not found in the file."
#define CONTOURXYZ_ZERO_AREA 468			// "X and Y values define a zero-area boundary: there is nothing to contour."
#define BAD_CONTOUR_LBL_FMT 469				// "Expected labelFormat=0, =1, =3, or =5."

#define XOP_EMPTY_CODE_FRAGMENT 470			// "The XOP's code fragment has nothing in it. It needs to be recompiled."
#define CANT_MODIFY_RGB_IMAGE 471			// "This is an RGB image; there is nothing to modify."
#define CANT_REPLACE_CONTOUR_TRACE 472		// "You can't replace a contour trace; it would just come back when the contours update."

#define OBJECT_DOES_NOT_EXIST 473			// "The named object does not exist."
#define WRONG_OBJECT_TYPE 474				// "The object is not of the specified type."

// End of Igor Pro 3.0 error codes.

#define FIRST_IGOR31_ERR 475				// Start of Igor Pro 3.1 error codes.

#define FUNCFIT_IND_VAR_MISMATCH	475		// "Number of independent variables in fit function does not match number of independent variable waves or columns."
#define FUNCFIT_DIMENSION_MISMATCH 476		// "Number of dimensions must match data wave."
#define FUNCFIT_ROW_MISMATCH 477			// "X waves must have the same number of rows as the data wave."
#define FUNCFIT_TOO_MANY_X_DIMS 478			// "Independent variable array must have no more than 2 dimensions."
#define FUNCFIT_ONE_X_MATRIX 479			// "When you use a matrix for independent variables, only one wave is allowed."
#define FUNCFITXYZ_XWAVE_REQUIRED 480		// "Your fitting function has more than one independent variable; you must specify at least one X wave using the /X flag."
#define FUNCFIT_TOO_MANY_IND_VARS 481		// "Too many independent variables; the limit is 16."
#define BAD_EWave 482						// "Epsilon wave does not match parameter wave (length or number type)"
#define FITFUNC_RETURNED_NAN 483			// "The fitting function returned NaN for at least one X value."
#define FITFUNC_RETURNED_INF 484			// "The fitting function returned INF for at least one X value."
#define MDFUNCFIT_TOO_MANY_IND_VARS 485		// "The fitting function has more than 4 independent variables."
#define MDFUNCFIT_TOO_FEW_IND_VARS 486		// "FitFuncMD requires at least 2 independent variables; use FitFunc instead."
#define MDFUNCFIT_IND_VAR_MISMATCH 487		// "The data wave must have one dimension for each independent variables."
#define BAD_RWave 488						// "Residual wave does not match parameter wave (length or number type)"
#define CONSTRAINT_HOLD_CONFLICT 489		// "Fitting parameter ^3 is being held and also constrained."
#define CONF_WRONG_WAVES	490				// "Wrong number of confidence band waves for the options selected."
#define CRVFIT_CONFLEVEL_NOT_IN_RANGE 491	// "Confidence level must be between 0 and 1, corresponding to 0 to 100% confidence levels."
#define CONSTRAINT_ILLEGAL 492				// "When parsing constraint expression:\n  \"^2\"\nreceived error message \"^3\""
#define CONSTRAINT_MULTIPLE_K 493			// "The constraint expression\n  \"^2\"\n has more than one fit parameter in a single term"  
#define CONSTRAINT_NONLINEAR 494			// "The constraint expression\n  \"^2\"\n lacks a conditional operator (< or >)."
#define CONSTRAINT_K_OUT_OF_RANGE 495		// "Fit parameter ^0 should be in range K0 to ^1 in constraint expression \"^2\"."
#define CONSTRAINT_NO_CONDITIONAL 496		// "Fit parameter ^1 is out of range in constraint expression\n  \"^2\"."
#define CONSTRAINT_ILLEGAL_OP_BEFORE 497	// "Illegal operator \'^0\' before fit parameter ^1 in constraint expression \"^2\"."
#define CONSTRAINT_ILLEGAL_OP_AFTER 498		// "Illegal operator \'^0\' after fit parameter ^1 in constraint expression \"^2\"."
#define MD_MISSING_XWAVE 499				// "The X wave was null or missing"
#define MD_MISSING_YWAVE 500				// "The Y wave was null or missing"
#define MD_MISSING_ZWAVE 501				// "The Z wave was null or missing"
#define MD_MISSING_TWAVE 502				// "The T wave was null or missing"
#define CURVEFIT_MISSING_XWAVE 503			// "X wave was null or missing"
#define CURVEFIT_MISSING_PWAVE 504			// "The parameter wave was null or missing"
#define CURVEFIT_MISSING_EWAVE 505			// "The Epsilon wave was null or missing"
#define CURVEFIT_MISSING_WWAVE 506			// "The weighting wave was null or missing"
#define CURVEFIT_MISSING_DWAVE 507			// "The destination wave was null or missing"
#define CURVEFIT_MISSING_RWAVE 508			// "The residual wave was null or missing"
#define CURVEFIT_MISSING_CMATRIX 509		// "The constraint matrix wave was null or missing"
#define CURVEFIT_MISSING_CVECTOR 510		// "The constraint vector wave was null or missing"
#define CURVEFIT_MISSING_CONFWAVE 511		// "A confidence band wave was null or missing"
#define CONSTRAINT_REQUIRES_TEXT_WAVE 512	// "Expected text wave containing constraint expressions"
#define CONSTRAINT_TEXT_WAVE_MISSING 513	// "Text wave for constraints was null or missing"

#define BAD_MISC_RECORD_VERSION 514			// "A miscellaneous data record is too new or corrupted."
#define CANT_CREATE_HOME_PATH 515			// "An error occurred while creating the home path."

#define BAD_PAGE_ORIENTATION 516			// "Expected 'Portrait' or 'Landscape'."
#define BAD_PRINT_PREFS	517					// "Page-Setup "

#define CANT_REMOVE_NORMAL_RULER 518		// "The Normal ruler can not be removed."

#define CANT_HANDLE_NON_NATIVE_PICT 519		// "This operation can't handle a <other platform> picture."
#define ERR_520 520							// Error 520 is available for duty.

#define SIN_FIFO_EXPECT_1CHAN 521			// "FIFO not setup properly (no channel info)"
#define SIN_FIFO_BAD_NUM_TYPE 522			// "FIFO number type not valid for sound input (8 and 16 bit ints only)"
#define NO_SIN_AVAIL 523					// "Sound input not available."
#define SIN_BAD_WAV_TYPE 524				// "Wave number type not valid for sound input (8 and 16 bit ints only)"
#define SIN_BAD_WAV_DIMS 525				// "Wave used for sound input can't have more than 2 dimensions."
#define SIN_FIFO_ALREADY_STARTED 526		// "Sound input to FIFO already started."
#define SIN_FIFO_ALREADY_STOPPED 527		// "Sound input to FIFO already stopped."
#define SIN_FIFO_ERROR  528					// "An error occurred during sound input to FIFO (probably overflow)."
#define TOOMANYCHANS 529					// "Sound input does not support the number of channels specified."
#define TOOMANYBITS 530						// "Sound input does not support the number of bits specified."
#define FREQNOTAVAIL 531					// "Sound input does not support the sampling rate specified."
#define SI_ERROR_SET_GAIN 532				// "Can't set the sound input gain."
#define SI_ERROR_SET_AGC 533				// "Can't set the sound input AGC."
#define AUDIO_FORMAT_NOT_AVAIL 534			// "The sound device does not support the specified sample rate, sample bits and/or channels."
#define AUDIO_BAD_VIBS 535					// "An audio related error occurred."

#define FAIL_READING_ENHMF 536				// "An error occurred while attempting to read an enhanced metafile."
#define FAIL_READING_AldusMF 537			// "An error occurred while attempting to read a placeable metafile."
#define FAIL_READING_WindowsMF 538			// "An error occurred while attempting to read a windows metafile."
#define FAIL_READING_DIB 539				// "An error occurred while attempting to read a device independent bitmap."

#define CANT_READ_THAT_GRAPHICS_TYPE 540	// "Can't read the graphics format of the specifed file."
#define NO_AUDIO_DEVICE 541					// "Could not find an audio device."
#define AUDIO_SYNCH_ONLY 542				// "Your audio setup does not support asynchronous output."
#define CLIP_NOT_AVAIL 543					// "Clipboard in use by another application."
#define PNG_WRITE_ERROR 544					// "An error occured while writing a PNG file."
#define CLIP_ERROR 545						// "A clipboard error occured."
#define EXPECT_SVAR 546						// "Expected name of string variable reference (SVAR)."
#define EXPECT_NVAR 547						// "Expected name of numeric variable reference (NVAR)."
#define EXPECT_0_OR_LPAREN 548				// "expected 0 or '('" */
#define EXPECT_KEYWORD 549					// "Expected keyword."

#define BAD_COLORTABLE_INDEX 550			// "IGOR color table index out of range"
#define BAD_COLORTABLE_HANDLE 551			// "Invalid IGOR color table handle"
#define BAD_COLORTABLE_PARAM 552			// "IGOR color table parameter out of range"

#define INCLUDE_FILE_ALREADY_OPEN_AS_NOTEBOOK 553	// "Included file is already open as a notebook"
#define INCLUDE_FILE_ALREADY_OPEN_AS_ANOTHER 555	// "Included file is already open as another type of window (e.g., help window)"

#define NOT_EXPERIMENT_FILE 555				// "This does not appear to be a valid IGOR experiment file"

#define FUNCLIST_UNKNOWN_OPTION 556			// "FunctionList does not recognize the option \"^0\""
#define CRVFIT_CONF_NO_MV 557				// "Confidence bands are not supported for multi-variate curve fits."
#define FIFO_ERR_SWAP_INFO 558				// "Probable bug while swapping FIFO info; please contact WaveMetrics tech support"
#define FIFO_ERR_SWAP_CHANINFO 559			// "Probable bug while swapping FIFO channel info; please contact WaveMetrics tech support"

#define ODE_STEP_SIZE_TOO_SMALL 560			// "The integration step size has gotten too small. You might try relaxing the error tolerance."
#define ODE_BAD_METHOD 561					// "The method keyword is set to a value that does not correspond to a supported method."
#define ODE_DERIVE_FUNC_MISSING 562			// "You must specify a function to calculate the derivatives."
#define ODE_PARAM_WAVE_MISSING 563			// "The pWave keyword is missing, or the wave specified by pWave keyword does not exist."
#define ODE_RESWAVE_MISSING 564				// "The resWave keyword is missing, or a wave or waves specified by the resWave keyword does not exist."
#define ODE_XWAVE_WRONG_NPNTS 565			// "The X wave (specified by the xvalues keyword) must have the same number of rows as the result waves (specified by the resWave keyword)."
#define ODE_SCALEWAVE_WRONG_NPNTS 566		// "The error scaling wave (errscale keyword) must have the same number of points as the parameter wave (pwave keyword)."
#define ODE_ZERO_SCALE 567					// "One or more of the values in the error scaling wave (errscale keyword) is zero."
#define ODE_MISSING_SCALE_WAVE 568			// "You have selected an error checking method (errorMethod keyword) that requires a wave to specify error scaling (errscale keyword)."
#define ODE_MISMATCHED_YWAVES 569			// "The lengths of the result waves (reswave keyword) must all be the same."
#define ODE_BAD_FUNCTION 570				// "The function ^0 has the wrong parameter types or the wrong number of parameters to use as an ODE function."
#define CURVEFIT_MISSING_YWAVE 571			// "Y wave was null or missing"
#define BAD_NOTEBOOK_PICTURE_MODE 572		// "Invalid notebook picture mode."
#define kHistEqHistWaveMustBe1D 573			// "The histogram wave must be a 1D wave"
#define kRequiresImageData 574				// "This operation works on image data (i.e., integer) waves only"
#define kImageHistSourceBadWaveType 575		// "Image histogram source wave must be 2D or 3D. If 3D, it must have exactly 3 layers."
#define kNumPointsMustBeEven	576			// "Histogram levels must be in pairs."
#define	kNoSuchThresholdMethod	577			// "No such threshold method."
#define kWavesMustBeSameSize	578			// "Both images must be the same size."
#define kWantsUnsignedCharData	579			// "This image operation supports only unsigned char (/B/U) data"
#define kMissingClipValue		580			// "Adaptive Histogram requires clip value (/C flag)"
#define kBadRegionSpecifier		581			// "Bad region specifier. Check that the image can be evenly divided into the number of specified regions."
#define kBadClipValue			582			// "Clipping value must be positive. /C=1 returns the original image."
#define kBadNumberOfBins		583			// "Bad value for the number of bins."
#define kNumHRegions2Big		584			// "The number of horizontal regions is too big."
#define kNumVRegions2Big		585			// "The number of vertical regions is too big."
#define kRequiresEqualRegions	586			// "Image size must be a multiple of the number of regions."
#define kRequiresMin4Regions	587			// "Adaptive histogram requires at least 4 sub regions."
#define kBadWaveInitialization	588			// "Bad wave initialization"
#define kBadWaveForMask			589			// "Bad wave for mask"
#define kNoSuchFilter			590			// "No such filter"
#define kNoSuchStructElement	591			// "Structure element is undefined."
#define kNeeds3DWave			592			// "This operation requires a 3D wave."	
#define kUnspecifiedThreshold	593			// "Threshold has not been specified"
#define kMethod3NotSupported	594			// "Method 3 is not supported for this operation."
#define kMustHaveROIWave		595			// "An ROI wave is required for this operation."
#define BAD_RECENT_FILES_PREFS	596			// "Recent Files "
#define INCLUDE_EXTENSION_NOT_ALLOWED 597	// "You must omit the \".ipf\" extension in #include statements."
#define PRINTER_DRIVER_SCREWUP 598			// "The printer driver returned an unexpected error. Check the driver and the printer."

#define FIRST_IGOR312_ERR 599				// Start of Igor Pro 3.12 error codes.
#define XOP_LOAD_LIBRARY_FAILED 599			// "Can not load XOP. It may be incorrectly constructed by the development system or may be corrupted."
#define CANT_FIND_XOP_MAIN 600				// "Can't find XOP's main function. (It must be declared HOST_IMPORT.)"
#define kBadROIDimensions		601			// "ROI dimensions do not match the target image."
#define BAD_FILE_TYPE_EXT_SPEC	602			// "A filetype/extension specifier must have exactly 4 characters or start with '.'"

#define FIRST_IGOR313_ERR 603				// Start of Igor Pro 3.13 error codes.

// These were created as platform-independent error codes that can be used when standard C file I/O routines return errors.
#define FILE_READ_ERROR 603					// "A file read error occurred."
#define FILE_WRITE_ERROR 604				// "A file write error occurred."
#define FILE_POS_ERROR 605					// "A file position error occurred."
#define FILE_OPEN_ERROR 606					// "The file can not be opened."
#define FILE_CLOSE_ERROR 607				// "An error occurred while closing the file."
#define FILE_PERM_ERROR 608					// "The file can not be opened because it is already open."
#define FILE_CREATE_ERROR 609				// "Error creating file. The file name or path may be invalid or the file may already exist."
#define FILE_EOF_ERROR 610					// "While reading a file, the end-of-file was reached."

#define XOP_LINK_FAILED 611					// "XOP dynamic link failed. The XOP may require a different version of Igor or may require additional DLLs."
#define SOME_LEVELS_MISSING 612				// "Some level crossings were not found."

// Root finder errors
#define ROOTS_TOO_MANY_EVALS 613			// "The root finder has performed more than the maximum allowed number of iterations."
#define ROOTS_NO_PROGRESS 614				// "The root finder not making progress. It may have gotten trapped at a point that is not a root, or your system may not have any roots."
#define ROOTS_NO_BRACKET 615				// "The root finder was unable to bracket a root before starting; for 1D roots it must find two X values where your function has opposite signs."
#define ROOTS_MISSING_X_WAVE 616			// "The X wave was not found."
#define ROOTS_MISSING_POLY_WAVE 617			// "The wave containing polynomial coefficients was not found."
#define ROOTS_MISSING_PWAVE 618				// "The parameter wave associated with the function ^0 was not found."
#define ROOTS_FUNCTION_TOO_MANY_PARAMS 619	// "Your function parameter wave had just one column, but your function has more than one independent variable. Your parameter wave must have a column for each independent variable."
#define ROOTS_FUNCTION_PARAMS_MISMATCH 620	// "The number of columns in the parameter wave must match the number of independent variables in the function."
#define ROOTS_WRONG_FUNCTION_PARAMS 621		// "You list ^0 functions; the functions must have ^1 parameters: a parameter wave and ^0 independent variables."
#define ROOTS_X_MISMATCH 622				// "Number of X values (length of X wave) does not match the number equations in your system."
#define ROOTS_FUNCTION_BAD_PARAMS 623		// "The function ^0 has the wrong parameter types or the wrong number of parameters."

#define CREATE_PROCESS_ERROR 624			// "CreateProcess failed with error \"^2\"."
#define CONTOUR_NODATA 625					// "Contour data is entirely blank (NaN). There is nothing to contour."

// End of Igor Pro 3.1 error codes.

#define IGORMENUMODE_BADMENUNAME 626		// "The menu \"^1\" is not recognized by SetIgorMenuMode."
#define IGORMENUMODE_BADITEMTEXT 627		// "The menu item text \"^1\" is not recognized by SetIgorMenuMode."
#define IGORMENUMODE_CANTDOTHATMENU 628		// "You can't enable/disable the items in that menu. Try disabling the item it is attached to in the parent menu."
#define IGORMENUMODE_BADACTION 629			// "The action to perform must be one of EnableItem, DisableItem, EnableAllItems or DisableAllItems."

#define CURVEFIT_BAD_POLY2D_ORDER 630		// "Poly 2D order must be at least 1."
#define CURVEFIT_MV_MISSING_XWAVES 631		// "You are fitting a multivariate fit function to a one-column data set.\nYou must select an X wave for each independent variable."
#define CURVEFIT_MISSING_CONF_LEVEL 632		// "You have checked the Error Analysis checkbox. You must enter a confidence level."
#define CURVEFIT_NO_Y_WAVE 633				// "You have not selected data to fit in the Y Data menu."
#define CURVEFIT_AMBIGUOUS_COEFS 634		// "Igor can't tell how many fit coefficients are required by the fit function. You must select a wave from the Coefficient Wave menu",
#define CURVEFIT_NEWWAVE_CONFLICT 635		// "You have selected 'New' from the %s menu, but a wave with the name %s already exists."
#define CURVEFIT_MISSING_NEWWAVE 636		// "You have selected 'New' from the %s menu, but you have not entered a name in the New Wave box."
#define CURVEFIT_MISSING_GUESS_USER 637		// "You have selected a user-defined fit function so you must enter an initial guess for every fit coefficient."
#define CURVEFIT_MISSING_GUESS_BUILTIN 638	// "You have selected manual guesses so you must enter an initial guess for every fit coefficient."
#define CURVEFIT_MISSING_GUESS_HOLD 639		// "You have checked a box to hold the value of a fit coefficient, but you have not entered a value for that coefficient. Go to the Coefficients Tab to do this."
#define CURVEFIT_MISSING_EPSILON_VALUE 640	// "You have selected an Epsilon wave so you must enter values in the Epsilon column of the Coefficients list."
#define CF_PLOTIT_IMSORRY 641				// "I'm sorry- it was impossible to plot your fitting function because"
#define CF_PLOTIT_MultiVariate 642			// "this service is not offered for multivariate fitting functions."

#define EXPECTED_GRAPH 643					// "expected name of a graph"

#define CF_PLOTIT_NoYWave 644				// "you have not selected a data wave from the 'Fit To' menu on the Input Data tab."
#define CF_PLOTIT_NoGraph 645				// "you need to make a graph."
#define CF_PLOTIT_WaveNotOnGraph 646		// "the fit data is not on the top graph."
#define CF_PLOTIT_YWaveIsXWave 647			// "on the top graph the fit data is used as an X wave."
#define CF_PLOTIT_ErrorMakingWave 648		// "there was an error making the destination wave. Igor is probably out of memory."
#define CF_PLOTIT_CouldntAddDisplay 649 	// "there was an error adding the destination wave to the graph. Igor is probably out of memory."
#define CF_PLOTIT_EquationNotAccepted 650	// "there was an error involving creating the function expression. Igor is probably out of memory."
#define CF_PLOTIT_NoCoefs 651				// "Igor couldn't determine the number if fit coefficients.\nTry selecting an appropriate coefficients wave from the Coefficients Wave menu."
#define CF_PLOTIT_ErrorMakingCWave 652 		// "Igor was unable to create a coefficients wave for the operation. Igor is probably out of memory."
#define CF_PLOTIT_NoGuess 653				// "you must enter initial guesses for every coefficient."
#define CURVEFIT_BADNCOEF 654				// "The fitting function requires a coefficient wave with ^0 points."
#define BAD_MWave 655						// "You have provided a mask wave that does not match your data wave."
#define CURVEFIT_MISSING_MWAVE 656			// "Mask wave was null or missing"

#define kMustHaveRowsAndCols	657			// "You must specify /N={extraRows,extraCols} for this operation" 	19MAR99		4.00D01
#define kBadCMAPWaveType		658			// "This wave is not appropriate as a CTAB."						19MAR99		4.00D01
#define kRowAndColsMustBeReasonable	659		// "The requested padding is not appropriate for this wave."		19MAR99		4.00D01
#define kBadWidthOrHeight			660		// "Improper width or height specification."
#define kBadXYWaves					661		// "Mismatch or inappropriate x, y waves."
#define kBadSeedPixel				662		// "Bad seed pixel specification."									09AUG99
#define kBadIPParamSpec				663		// "Bad IP Parameter specification."								09AUG99
#define kExpectedImageTypeName		664		// "Expected name of image file type."
#define kBadSrcWave					665		// "Source wave is bad."								AG	07JUN01
#define kXYZWaveMustBeSP			666		// "XYZ waves must be of type SP (NT_FP32)"				AG	28JUL99
#define kBadMultipleImageCount		667		// "Bad count for multiple images."						AG 	20SEP99
#define kBadValueForFirstImage		668		// "Bad value for first image."							AG	20SEP99
#define kWantsNewerQuickTime		669		// "Operation requires a newer version of QuickTime."	AG 	20SEP99

#define kMustSpecifyDataWave		671		// "Data wave must be specified (see /D flag)."			AG 	11OCT99 
#define kIncompatibleFlagOptions	672		// "Flag options appear to be incompatible."			AG 	29OCT99
#define kAGBloatedFortranCrap		673		// "AG's bloated fortran crap is not compiled."			AG	01NOV99
#define kExpectSquareMatrix			674		// "Expected square matrix."							AG 	05NOV99
#define kInsufficientPointsInWave	675		// "Insufficient number of points in the input wave."	AG	12NOV99
#define kAllPointsColinear			676		// "All points are co-linear."							AG 23NOV99
#define kAllPointsCoPlanar			677		// "All points are co-planar."							AG 23NOV99
#define kFailedConsistency			678		// "Failed consistency test."							AG 23NOV99
#define kNotConvex					679		// "Data does not represent convex set."				AG 23NOV99
#define kFailedEulerTest			680		// "Data failed Euler test."							AG 23NOV99
#define kExpectedTripletWave		681		// "Expected 3 column (Triplet) wave."					AG 03MAY00
#define kExpectedNT_FP32			682		// "Expected an SP (single precision float) wave."		AG 24MAY00
#define kExpectedSPorDP				683		// "Expected SP or DP wave"
#define kInputTooLarge				684		// "Input is too large"
#define kExpectNoNaNorInf			685		// "Source wave should not contain NAN's or INF's"		AG 18MAY01
#define kFailedInternalConsistencyTest	686	// "Failed internal consistency text"

#define ButtonRecursionDetected	687		// "Recursion prevented. The expression \"^0\" causes itself to be executed again, and needs to be changed."
#define NO_SUCH_TOOL_NAME	688			// "Expected \"normal\", \"arrow\", \"text\", \"line\", \"rect\", \"rrect\", \"oval\", or \"poly\"."

#define NO_MOVIE 689					// "no movie"
#define FAILED_TO_PLAY_MOVIE 690		// "failed to play movie"
#define EXPECT_WATERFALL 691			// "expected a waterfall plot"
#define X_VECTOR_MISMATCH 692			// "x vector mismatch"
#define Y_VECTOR_MISMATCH 693			// "y vector mismatch"
#define EXPECTED_INSTANCE 694			// "expected instance"
#define UNMATCHED_CONDITIONAL 695		// "unmatched ?: conditional"
#define NOT_IN_MACROS 696				// "this syntax is not allowed in macros -- use functions"
#define LINK_NO_XFUNC 697				// "During link, couldn't find external function."
#define LINK_TYPE_XFUNCMISMATCH 698		// "During link, external function did not match."
#define LINK_NO_FUNC 699				// "During link, couldn't find user function."
#define LINK_TYPE_MISMATCH 700			// "During link, user function did not match."
#define LINK_NO_CONST 701				// "During link, couldn't find constant."
#define SEMICOLON_EXPECTED 702			// "Expected semicolon."
#define EXPECT_OBJ_NAME 703				// "Expected object name."
#define EXPECT_CONSTANT_OR_LITERAL 704	// "Expected symbolic constant or literal"
#define EXPECT_FUNC_NAME 705			// "Expected function name."
#define FUNCREF_TYPE_INCONSISTENT 706	// "Function reference type inconsistent."
#define CANT_USE_FUNCREF_HERE 707		// "Can't use a function reference here."
#define EXPECT_LOCALVAR_NAME 708		// "Expected a local variable name."
#define REF_VAR_DIFFERENT_TYPE 709		// "Reference variable is of a different type."
#define COULD_NOT_FIND_PROTO_FUNC 710	// "Couln't find prototype function."
#define EXPECT_FUNC_REF 711				// "Expected function reference."
#define NO_STATIC_FUNC_HERE 712			// "Can't use a static function here."
#define NO_PROMPT_THIS_TYPE 713			// "Can't prompt for this type of variable."
#define EXPECT_POPUP 714				// "Expected popup keyword."
#define NO_PROMPT_DEFINED 715			// "No prompt defined for this variable."
#define FASTOP_TOO_MANY_PRODS 716		// "Too many product terms."
#define FASTOP_SYNERR 717				// "Syntax does not conform to FastOp requirements."
#define WAVE_LENGTH_MISMATCH 718		// "Wave lenght mismatch."
#define WAVE_TYPE_MISMATCH 719			// "Wave type mismatch."
#define COMPLEX_INT_NOT_AVAIL 720		// "Complex integers not supported here."
#define COMPLEX_TO_REAL_LOSS 721		// "Complex wave used in real expression."
#define DUP_CONST 722					// "Duplicate constant."
#define WRONG_IGOR_VERS 723				// "A more recent version of Igor is required."
#define NOGRAF_OR_PANEL_OR_TABLE 724	// "Expected graph, panel or table."
#define EXPECT_WIN_NAME 725				// "Expected window name."

#define NOSUBRANGECATWAVE 726				// "Subranged category waves is not currently supported."
#define ONLY_ONE_RANGE 727					// "Only one dimension may have a range."
#define DUPLICATE_RESOURCES 728				// "This experiment can not be loaded in Carbon Igor because of duplicate resources. Open and save it in pre-Carbon Igor to fix the problem."
#define EXPECT_GUDE_NAME 729				// "Expected guide name."
#define NO_REDEFINE_BI_GUIDE 730			// "Can't redefine a built-in guide."
#define NO_SUCH_GUIDE 731					// "Specified guide does not exist."
#define NO_MIX_GUIDE_ORIENTATION 732		// "Guide orientation mismatch."
#define NO_SWITCH_GUIDE_ORIENTATION 733		// "Can't switch guide orientation."
#define ILLEGAL_DUAL_GUIDE_POS 734			// "Illegal dual guide position."
#define NO_CONTROLS_IN_SUBGRAPH 735			// "Can't put controls in a subgraph."
#define NO_PANEL_IN_SUBGRAPH 736			// "Can't put panels in a subgraph."
#define WRONG_GUIDE_ORIENTATION 737			// "Wrong guide orientation."
#define GUIDE_IN_GRAPH_ONLY 738				// "Guide is for graphs only."
#define NOT_VALID_GUIDE_NAME 739			// "Invalid guide name."
#define NO_SUCH_SUBWINDOW 740				// "Specified subwindow does not exist."
#define NOT_AVAILABLE_ON_SUBRANGE 741		// "Action not available on subrange."
#define COMPILE_ONLY 742					// "This feature can only be used in user functions."

#define MACROLIST_UNKNOWN_OPTION 743		// "MacroList does not recognize the option \"^0\""
#define EXPECTED_GRAF_OR_PANEL 744			// "Expected graph or panel name."

#define SEARCH_CANT_FIND_WINDOW	745			// "Can't find the referenced window in the experiment. Try doing the search again."
#define SEARCH_TOO_MANY_TERMS 746			// "Too many search terms - only 8 are allowed."
#define SEARCH_FILE_MODIFIED 747			// "The file was modified after the search. Try doing the search again."

#define WAVE_REF_FAILED	748					// "Failed to resolve a local WAVE \"^0\" reference to a global wave."
#define NVAR_REF_FAILED 749					// "Failed to resolve a local NVAR \"^0\" reference to a global variable."
#define SVAR_REF_FAILED 750					// "Failed to resolve a local SVAR \"^0\" reference to a global string."

#define DIM_LABEL_TOO_LONG 751				// "Dimension labels are limited to 31 characters."

#define CF_PLOTIT_NoGuessBuiltin 752		// "you must enter function coefficient values. Select Manual from the Guess Method menu first."
#define CF_PLOTIT_NoGuessPolyLine 753		// "you must enter function coefficient values. For poly and line fits, check the Hold box in order to enter values."

#define COLUMN_INFO_UNKNOWN_SPECIFIER 754	// "Expected 'C', 'F', 'W', or 'N' in column info specifier."
#define COLUMN_INFO_EXPECTED_NAME 755		// "Expected a wave name in column info specifier."
#define COLUMN_INFO_EXPECTED_NUMBER 756		// "Expected a number in column info specifier."
#define COLUMN_INFO_BAD_NAME 757			// "A name in the column info specifier contained illegal characters."
#define COLUMN_INFO_BAD_NAME_TERM 758		// "Missing comma or semicolon after a name in the column info specifier."
#define COLUMN_INFO_NAME_TOO_LONG 759		// "A name in the column info specifier can not exceed 31 characters in length."
#define COLUMN_INFO_BAD_NUMBER 760			// "Bad number in the column info specifier."
#define COLUMN_INFO_BAD_NUMBER_TERM 761		// "Missing comma or semicolon after a number in the column info specifier."
#define BAD_FIXED_FIELD_NUMBER_OF_COLUMNS 762	// "The number of columns in a fixed field file must be between 1 and 10000."
#define BAD_FIXED_FIELD_FIELD_WIDTH 763			// "The field width in a fixed field file must be >= 1."
#define EXPECT_IMAGE_CONTOUR_TRACE_NAME 764		// "Expected name of image, contour, or trace in top graph."
#define CL_LOOKUP_REQUIRES_CTAB 765				// "ColorScale lookup requires ctab keyword."
#define CURVEFIT_NOSTATICFUNCTIONS 766			// "Static function references are not allowed as curve fitting functions."
#define EXPECTED_GRAPH_TABLE_LAYOUT_PROCEDURE_NOTEBOOK_CMDWIN_BUT_NOT_PANEL 767	 // "This operation is for graphs, tables, layouts, notebooks, procedure windows, or the command/history window."
#define EXPECTED_TARG_PROC_CMDWIN_NAME 768		// "expected name of a target window, procedure window, or \"kwCmdHist\"."
#define CURVEFIT_NOTENOUGHPOINTS 769			// "You must have at least as many data points as fit parameters."
#define BAD_TIMEUNIT 770						// "Expected the name of a time unit like sec, week, year, etc."
#define MANDATE_INCMUSTBEINTEGER 771			// "The manual tick increment for a date/time axis must be an integer."
#define DATEAXIS_NOSWAPENDS 772					// "Date/time axes do not support reversed scaling."
#define FLAG_MUST_BE_FIRST_OR_AFTER_SLASH_W 773	// "this must be first flag in command, or immediately follow /W."
#define EXPECTED_TARG_PROC_NAME 774				// "expected name of a target window or procedure window."

#define EXPECTED_GRAPH_MARGIN_VAL 775			// "Expected a graph margin value."
#define EXPECTED_GRAPH_PLOT_AREA_VAL 776		// "Expected a graph plot area width or height value."

#define LAYOUT_CAN_NOT_BE_COMPILED 777			// "The Layout operation can not be used in a function. Use NewLayout instead."
#define APPENDTOLAYOUT_CAN_NOT_BE_COMPILED 778	// "The AppendToLayout operation can not be used in a function. Use AppendLayoutObject instead."
#define BAD_FRAME_VALUE 779						// "Expected a frame value: 0=none, 1=single, 2=double, 3=triple, 4=shadow."
#define BAD_LAYOUT_OBJECT_TYPE 780				// "Expected a page layout object type keyword: graph, table, picture or textbox."
#define EXPECTED_LAYOUT_NAME 781				// "Expected the name of a page layout window."
#define LAYOUT_USE_TEXTBOX_CMD 782				// "Can't append a textbox via AppendLayoutObject. Use Textbox or Legend instead."
#define FLAG_ALLOWED_JUST_ONCE 783				// "This flag can be used only once per command."

#define OPTIMIZE_NOXFLAG 784					// "When optimizing a multivariate function, you must provide a starting guess using the /X flag."
#define OPTIMIZE_NAN 785						// "The user function you are trying to optimize has returned NaN (Not a Number)."
#define OPTIMIZE_NOBRACKET 786					// "The Optimize operation was unable to find a pair of X values that bracket the minimum (or maximum). Use /L and /H to specify bracketing values."
#define OPTIMIZE_NOPROGRESS 787					// "The Optimize operation could not find a better solution than your starting guess. Your function may be too non-linear, the stopping tolerance may be too large, or your starting guess is a solution."
#define OPTIMIZE_TOOMANYITERATIONS 788			// "The Optimize operation has performed more ^0 iterations, which is more than the maximum allowed."
#define OPTIMIZE_MAXSTEPSIZE 789				// "The Optimize operation  has exceded the maximum step size. It may be that your function is unbounded or approaches a value assymptotically."
#define OPTIMIZE_NTYPSIZEMISMATCH 790			// "The number of values used with the /R flag must match the number of X values."
#define OPTIMIZE_CRITICALPOINT 791				// "Your starting guess is too near a critical point and Optimize can't procede. Try a different starting guess."

#define EXPECTED_LOCAL_NUM_OR_STR_VAR_NAME 792	// "Expected the name of a local numeric variable or NVAR, or a local string variable or SVAR."
#define SSCANF_RAN_OUT_OF_CONVERSIONS 793		// "The sscanf format string does not have enough conversion specifiers or there are too many output variables."
#define SSCANF_TOO_MANY_CONVERSIONS 794			// "The sscanf format string has too many conversion specifiers or there are too few output variables."
#define SSCANF_REQUIRES_L 795					// "\"%e\", \"%f\", and \"%g\" are not allowed in sscanf. Use "\"%le\", \"%lf\", and \"%lg\" instead."
#define SSCANF_L_NOT_ALLOWED 796				// "Do not put an 'l' after the '%' in an sscanf format string."
#define SSCANF_UNKNOWN_FORMAT 797				// "sscanf unsupported format."
#define SSCANF_IN_FUNCTIONS_ONLY 798			// "sscanf can be used in user functions only."
#define SSCANF_EXPECTED_NUM_GOT_STR 799			// "sscanf expected a numeric variable and got a string variable. Check the format string and the variable list."
#define SSCANF_EXPECTED_STR_GOT_NUM 800			// "sscanf expected a string variable and got a numeric variable. Check the format string and the variable list."
#define SSCANF_SCANSET_NOT_TERMINATED 801		// "A scanset (\"%[...]\") was used in sscanf but there was no trailing \"]\" character."
#define SSCANF_TOO_MANY_PARAMS 802				// "sscanf is limited to 100 localVar parameters."
#define CONSTRAINTS_NOT_AVAILABLE 803			// "This version of Igor does not support curve fitting with constraints."
#define EXPECTED_NUMERIC_WAVE 804				// "Expected numeric wave."
#define EXPECTED_TEXT_WAVE 805					// "Expected text wave."

#define CURVEFIT_BADAUTODESTLEN 806				// "The Destination length must be 2 or greater, or \"Auto\" or zero."
#define CF_PLOTIT_NoAllAtOnce 807				// "this service is not offered for all-at-once fitting functions."

// End of Igor Pro 4-compatible error codes.

#define FIRST_IGOR5_ERR 808

#define CANT_OPEN_SPECIFIED_PRINTER_DRIVER 808	// "Unable to open the printer driver for this page setup record."
#define INCOMPATIBLE_PRINT_RECORD 809			// "The page setup record is not compatible with the current printer driver."

#define GETPATH_CALLBACK_OBSOLETE 810			// "The GetPath XOP callback is no longer supported. The XOP needs to be updated for this version of Igor."
#define TU_CALLBACKS_OBSOLETE 811				// "The TUWriteFile and TUInsertFile XOP callbacks are no longer supported. The XOP needs to be updated for this version of Igor."

#define NOT_IMPLEMENTED 812						// "This feature is not yet implemented."
#define FCMD_BAD_INTERACTIVE 813				// "In /I=i, expected i between 0 and 3."
#define FCMD_CANT_COPY_SELF 814					// "Can't copy a file onto itself."
#define FCMD_CANT_MOVE_SELF 815					// "Can't move a file onto itself."

#define XOP_68K_NOT_SUPPORTED 816				// "This is a 68K XOP. It is not compatible with this version of Igor."

#define EXPECTED_MENU_ITEM_FLAG 817				// "Expected a menu item flag (/Q)."

#define XOP_CANT_RUN_ON_OS_X 818				// "The XOP '^3' can not run on Mac OS X."

#define NEGATIVE_FIFO_SIZE 819					// "The fifo size is negative, which is not allowed. Perhaps you used a null variable to set the size."
#define FIFO_SIZE_TOO_BIG 820					// "The fifo size is too big."

#define CANT_TARGET_ROOT_DIRECTORY 821			// "This operation does not permit targeting the root directory of a volume. Target a sub-directory."

#define kBadParameters	822						// "One or more parameters are inappropriate."
#define kParameterOutOfRange	823				// "Parameter is out of range"
#define kBadStructureElement	824				// "Bad Structure Element"
#define kExpectedComplexWave	825				// "Expected Complex Wave"
#define kDivideByZero			826				// "Operation failed beacue of a divide by zero condition"
#define kBadParam				827				// "Bad parameter."
#define kBadWaveletSpecification	828			// "Bad Wavelet Specification."
#define kBadOutputSpecification		829			// "Bad Output Specification."
#define kBadOffsetsSpecification	830			// "Bad Offset Specification."
#define kBadNumCoeff				831			// "Bad Number of Coefficients."
#define kBadRangeSpecification		832			// "Bad range specification."
#define kDestinationMustBeSpecified 833			// "Destination wave must be specified."
#define kRequirePositiveParameter   834			// "Positive number required."
#define kBadTriangulationWave		835			// "Bad Triangulation Wave."
#define kBadDestinationWave			836			// "Bad Destination Wave."
#define kDoesNotSupport4D			837			// "Does not support 4D waves."
#define kBadUserFunctionFormat		838			// "Bad user function format."
#define FFT_COLS_EVEN				839			// "The number of columns must be even."
#define kWave_Scaling_Mismatch		840			// "Wave Scaling Mismatch"
#define kBadMatrixOPToken			841			// "Bad MatrixOPs token."
#define kMatrixDimMismatch			842			// "Matrix dimensions mismatch."
#define kNeedSquareMatrix			843			// "Expected square matrix."
#define kRequireLeftHand			844			// "Left Hand Side is required."
#define kLeftSideMustBeWaveOrVar	845			// "Left Hand Side must be a wave."
#define kBadTokenIndex				846			// "Bad matrix token index."
#define kMatrixStackOutOfSync		847			// "Matrix stack is out of sync."
#define kUnbalancedParenthesis		848			// "Unbalanced parenthesis"
#define kBadProcessRange			849			// "Bad process range."
#define kExpectedOperator			850			// "Expected operator."
#define kCouldNotFindData			851			// "Could not find data."
#define kNaNsNotAllowed				852			// "NaNs are not allowed in this operation."
#define kExpectedPrefixOperator		853			// "Expected prefix operator."
#define kUnknownPrefixOperator		854			// "Unknown prefix operator."
#define kBadDataType				855			// "Bad data type."
#define kExpectRealMatrix			856			// "Expected real matrix."
#define kCantSubtractMatrixFromScalar 	857		// "Can't subtract matrix from a scalar."
#define kBadROIWaveType					858		// "Bad ROI wave type."
#define kBadNumberOfClasses				859		// "Bad number of classes."
#define kLAPACKError					860		// "LAPACK returned an error.  Check your input matrix."
#define kNoComplexSupport				861		// "Does not support complex waves or variables."
#define kAPMMemoryError					862		// "Arbitrary Precision Math memory error."
#define kDoesNotSupportNaNorINF			863		// "Does not support NaN or INF."
#define kGeneralException				864		// "IGOR encountered an insufficient memory condition or an unknown system exception."
#define kBadSnakeWaves					865		// "Bad snake wave specification."
#define kInsufficientSnakePoints		866		// "Snakes must have a minimum of 4 points."
#define kBadGradImageWave				867		// "Bad gradient image wave."
#define kBadExternalEnergyWave			868		// "Bad external energy wave."
#define kBadNumberOfSamples				869		// "Bad number of samples."
#define kExpectedNT_FP64				870		// "Expected an DP (double precision float) wave."		

#define SA_NEED_3_COLUMN_WAVE	871				// "A 3-column wave is required."
#define SA_NEED_2_COLUMN_WAVE	872				// "A 2-column wave is required."
#define SA_MISSING_STEP_SIZE_WAVE 873			// "The simulated annealing step size wave is null or missing."
#define SA_MISSING_MINMAXXWAVE 874				// "The simulated annealing X limit wave is null or missing."
#define PROGRAMMED_STOP 875						// "Termination requested by function."
#define MISMATCHED_NUMBER_OF_POINTS 876			// "Different number of points for X and Y waves." used by New Graph dialog, DisplayWaves.cpp
#define EXPECT_1D_WAVE_FROM_LIST 877			// "2D waves must have a single row or column selected. Click the Add button and select the row or column in the list below."

#define OH_UNKNOWN_PARAM_TYPE 878				// "Valid parameter types are number, string, name, wave."
#define OH_EXPECTED_POSTFIX_CHAR 879			// "Expected trailing ), ] or }."
#define OH_TOO_MANY_PARAMS 880					// "The template contains too many parameters."
#define OH_RUNTIME_STRUCT_TOO_SMALL 881			// "The template requires a runtime structure larger than the specified maximum."

#define OH_TOO_MANY_FLAG_PARAM_GROUPS 882		// "The template contains too many flag parameter groups."
#define OH_TOO_MANY_MAIN_PARAM_GROUPS 883		// "The template contains too many main parameter groups."
#define OH_OPTIONAL_PARAMS_MUST_BE_AT_END 884	// "Optional flag or keyword parameters must appear at the end of a group introduced by (, [ or {."
#define OH_OPTIONAL_SIMPLE_MAIN_PARAMS_MUST_BE_AT_END 885	// "Optional simple main parameters must be the last parameter group in the template."
#define OH_BAD_NUMBER_OF_OPTIONAL_PARAMS 886	// "Number of optional parameters must be between 1 and 100."
#define OH_EXPECTED_EQUALS 887					// "Expected '='."
#define OH_EXPECTED_COMMA 888					// "Expected comma between parameters."
#define OH_EXPECTED_LEFT_PAREN 889				// "Expected '('."
#define OH_EXPECTED_RIGHT_PAREN 890				// "Expected ')'."
#define OH_EXPECTED_LEFT_BRACKET 891			// "Expected ']'."
#define OH_EXPECTED_RIGHT_BRACKET 892			// "Expected '['."
#define OH_EXPECTED_LEFT_BRACE 893				// "Expected '{'."
#define OH_EXPECTED_RIGHT_BRACE 894				// "Expected '}'."
#define OH_CANT_MIX_SIMPLE_AND_KEYWORD_PARAMS 895	// "You can't mix simple and keyword main parameters in the same operation."
#define OH_BAD_XOP_OPERATION_NAME 896			// "The XOP does not add an operation with this name."
#define OH_OPERATION_NOT_REGISTERED 897			// "The operation is not registered with Igor."
#define OH_COMPILABLE_BIT_MUST_BE_SET 898		// "The operation's compilable bit must be set in the XOPC resource."
#define OH_BAD_RUNTIME_PARAM_STRUCT_SIZE 899	// "The runtime parameter structure size is too small."
#define OH_OUT_OF_MEMORY 900					// "Operation Handler ran out of memory while parsing command template."
#define OH_SYNTAX_REQUIRES_PREFIX_CHAR 901		// "This type of optional parameter syntax works only with {}, [] or () syntax."
#define OH_TOO_MANY_LEFT_BRACKETS 902			// "Use only one pair of brackets to indicate optional parameters within {}, [] or () syntax."
#define OH_EXPECTED_NUMBER 903					// "Expected a number or numeric expression"
#define OH_EXPECTED_STRING 904					// "Expected a string or string expression"
#define OH_EXPECTED_NAME 905					// "Expected a name"
#define OH_EXPECTED_VARNAME 906					// "Expected name of a variable, NVAR or SVAR"
#define OH_EXPECTED_WAVE 907					// "Expected a wave"
#define OH_EXPECTED_WAVE_TYPE 908				// "Expected wave type (real, complex, text)"

#define NO_SUCH_OPERATION 909					// "There is no operation with this name."
#define BAD_OPERATION_VAR_NAME 910				// "Names of variables created by operations must start with V_ (numeric) or S_ (string)."

#define EXPECTED_NUM_WAVE 911					// "Expected a numeric wave."
#define EXPECTED_NUM_VAR_OR_NVAR 912			// "Expected the name of a numeric variable or an NVAR."

// These are used for XOP's calling user functions.
#define BAD_COMPILATION_INDEX 913				// "A call to an internal Igor routine (CallUserFunction) was made using stale information."
#define REQUIRES_NUMERIC_PARAMETER 914
#define REQUIRES_PASS_BY_REFERENCE_NUMERIC_PARAMETER 915
#define REQUIRES_COMPLEX_NUMERIC_PARAMETER 916
#define REQUIRES_PASS_BY_REFERENCE_COMPLEX_NUMERIC_PARAMETER 917
#define REQUIRES_STRING_PARAMETER 918
#define REQUIRES_PASS_BY_REFERENCE_STRING_PARAMETER 919
#define REQUIRES_NUMERIC_WAVE_PARAMETER 920
#define REQUIRES_COMPLEX_NUMERIC_WAVE_PARAMETER 921
#define REQUIRES_TEXT_WAVE_PARAMETER 922
#define REQUIRES_FUNCTION_REFERENCE_PARAMETER 923
#define UNKNOWN_PARAMETER_TYPE 924
#define REQUIRES_NUMERIC_RETURN_TYPE 925
#define REQUIRES_COMPLEX_NUMERIC_RETURN_TYPE 926
#define REQUIRES_STRING_RETURN_TYPE 927
#define UNKNOWN_RETURN_TYPE 928
#define FUNCTION_HAS_TOO_FEW_PARAMETERS 929
#define FUNCTION_HAS_TOO_MANY_PARAMETERS 930
#define PASS_BY_REF_BAD_PARAM_TYPE 931			// "Only numeric and string parameters can be pass-by-reference"

#define XOP_MACH_CANT_FIND_EXECUTABLE 932		// "Can't find the executable file in the XOP package."
#define XOP_MACH_CANT_LOAD_EXECUTABLE 933		// "Can't load the executable file in the XOP package."
#define XOP_MACH_CANT_FIND_MAIN 934				// "Can't find the main function in the XOP package."
#define XOP_MACH_CANT_RUN_ON_OS9 935			// "Mach XOPs run on OS X only, they can not execute on Mac OS 9."

#define STOP_RECURSING_THROUGH_FOLDERS 936		// This is not really an error. It is a signal to RecurseThroughFolders.

#define EXPECTED_WSIZE 937

#define WARN_OVERRIDE_NOT_IN_MAIN_PROC 938		// "Override functions should be in the main Procedure window to work reliably."
#define WARN_FUNCTION_IS_OVERRIDDEN 939			// "This function will be overridden by a pre-existing Override function."
#define STATIC_NEEDS_MODULE_NAME 940			// "This procedure window lacks the #pragma Module statement static functions require."

#define EXPECT_STRUCT_NAME 941					// "Expected structure name."
#define EXPECT_STRUCT 942						// "Expected structure."
#define EXPECT_COMPAT_STRUCT 943				// "Expected compatible structure."
#define EXPECT_WAVE_FIELD 944					// "Expected WAVE field."
#define EXPECT_NVAR_FIELD 945					// "Expected NVAR field."
#define EXPECT_SVAR_FIELD 946					// "Expected SVAR field."
#define EXPECT_FUNCREF_FIELD 947				// "Expected FUNCREF field."
#define EXPECT_LOCALCONSTANT_OR_LITERAL_EXPR 948	// "Expected locally defined constant or literal expression."
#define NOT_OPTPARAM_NAME 949					// "Expected name of optional parameter."
#define NO_OPTPARAM_FUNCS 950					// "Functions with optional parameters not allowed here."
#define DUP_PROCPICT 951						// "Duplicate procedure Picture."
#define EXPECTED_END 952						// "Expected End."
#define EXPECT_ASCII85Begin 953					// "Expected ASCII85Begin."
#define DUP_PROCSTRUCT 954						// "Duplicate procedure Structure."
#define STRUCT_FIELD_ARRAY_OUTOFBOUNDS 955		// "Illegal structure field array size."
#define NO_OVERRIDE 956							// "Can't use Override here."
#define NO_SUCH_STRUCT 957						// "No such structure exists."
#define STRUCTS_TOO_DEEP 958					// "Structure nesting too deep."
#define STRUCT_FIELD_INDEX_OB 959				// "Structure field index out of bounds."
#define EXPECT_STRUCT_FIELD 960					// "Expected structure field."
#define ILLEGAL_FIELD_FOR_FBIN 961				// "Illegal field in structure (String, NVAR etc.) for this use."
#define DUP_FIELD_NAME 962						// "Duplicate field name."
#define CANT_CHANGE_LOCKED_WAVE 963				// "Can't change a locked wave."
#define AXIS_NAME_USED 964						// "An axis of that name already exists."
#define BAD_MASTER_AXIS 965						// "Specified master axis not found."
#define BAD_AXIS_HOOK 966						// "Axis hook function not found."
#define WRONG_FUNC_PARAMS 967					// "Function input parameters or output not valid for this use."
#define AXIS_IN_USE 968							// "Axis is in use."
#define STRUCT_REF_ONLY 969						// "Structure input parameters must be pass-by-reference ('&' needed.)"
#define USING_NULL_REFVAR 970					// "attempt to use uninitialized pass-by-reference variable"
#define FONT_ERR 971							// "General font error."
#define FONT_CURVE_EXTRACT_ERR 972				// "Error extracting font outline curves."
#define INCOMPATIBLE_STRUCT_VERSION 973			// "Incompatible structure. The calling function is too old or too new for the called function."
#define STRUCT_ONLY_IN_FUNCTION 974				// "Structure parameters can be used only in user-defined functions"
#define BAD_FUNCREF 975							// "FUNCREF does not reference valid function."
#define BAD_WAVE_LIST 976						// "There was a problem in a list of waves."
#define REQUIRES_STRUCTURE_PARAMETER 977		// "Requires structure parameter". Used for XOP's calling user functions.
#define NOT_IN_THREADSAFE 978					// "Not allowed in ThreadSafe functions."
#define NOT_YET_IN_THREADSAFE 979				// "Not yet available in ThreadSafe functions."
#define INVALID_THREAD_GROUP 980				// "Invalid Thread Group ID or index."
#define WAVE_USED_BY_THREAD 981					// "Wave is in use by preemtive thread. Can't be resized or killed."
#define ILLEGAL_THREAD_PARAM 982				// "Parameter not allowd when spawning a preemptive thread."
#define ONLY_THREADSAFE 983						// "Function must be ThreadSafe."
#define NO_CALLS_OUTSIDE_IM 984					// "Functions in Independent Module can't call outside."

#define OH_BAD_STRUCT_SIZE 985					// "The structure is the wrong size."
#define OH_BAD_STRUCT_TYPE_NAME 986				// "This is the wrong type of structure."
#define OH_BAD_STRUCTURE_TYPE 987				// "Expected structure type, 0 or 1."
#define OH_EXPECTED_LITERAL_INTEGER 988			// "Expected a literal integer number."
#define FILE_STRUCT_MISMATCH 989				// "The file size does not match the structure size."

#define EXPECT_GREATER_THAN 990					// "Expected number greater than ^2"
#define FCMD_AS_NOT_ALLOWED 991					// "This command doesn't accept an \"as\" parameter."

#define NO_COMPLEX_TEXT_WAVES 992				// "Igor does not support complex text waves."

#define STRING_ACCESS_ON_NUMERIC_VARIABLE 993	// "An attempt was made to treat a numeric variable as if it were a string variable."
#define NUMERIC_ACCESS_ON_STRING_VARIABLE 994	// "An attempt was made to treat a string variable as if it were a numeric variable."
#define BAD_VARIABLE_DATA_TYPE 995				// "A variable data type must double, double complex, or string."
#define COMPILE_FAILED 996						// "Procedure compilation failed."
#define BAD_WIN_TYPE 997						// "This operation does not apply to this type of window."
#define BAD_PAGESETUP_SCALE 998					// "Expected a page setup scaling value in percent between 5 and 5000."
#define FLAG_ALLOWED_ONLY_ONCE 999				// "This operation allows each flag to be used only once in a single command."
#define KEYWORD_ALLOWED_ONLY_ONCE 1000			// "This operation allows each keyword to be used only once in a single command."

#define	BAD_NUM 1001							// "expected number"
#define	NAM_SYMB_BAD 1002						// "unknown/inappropriate name or symbol"
#define	LPAREN_MISSING 1003						// "expected left parenthesis"
#define	RPAREN_MISSING 1004						// "expected right parenthesis"
#define	OPAND_MISMATCH 1005						// "expected operand"
#define	OPTOR_MISMATCH 1006						// "expected operator"
#define	OPTOR_OPAND_MISMATCH 1007				// "operator/operand mismatch"
#define WRONG_NO_PARAMS 1008					// "wrong number of parameters"
#define BI_NO_LPAREN 1009						// "expected left parenthesis"
#define RPAREN_EXPECTED 1010					// "expected right parenthesis"
#define NON_NULL_PARAMS 1011					// "this function takes no parameters"
#define BI_NO_FCTN_THIS_NT 1012					// "function not available for this number type"
#define AMBIGUOUS_WAVE_POINT 1013				// "ambiguous wave point number"
#define BI_BAD_WAVEFORM 1014					// "expected wave name"
#define NO_SUCH_CURSOR 1015						// "cursor ^0 is not on the graph \"^1\""
#define NAM_SYMB_BAD_NO_CSR 1016				// "expected cursor name (A or B)"
#define LINE_TOO_LONG 1017						// "line too long"
#define RBRAKET_EXPECTED 1018					// "expected ']'"
#define LBRAKET_EXPECTED 1019					// expected '['
#define NO_USER_VARS_IN_FUNC 1020				// user variables are not allowed in functions

// Automation errors
#define AUTOMATION_RESERVED_1021 1021			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1022 1022			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1023 1023			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1024 1024			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1025 1025			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1026 1026			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1027 1027			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1028 1028			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1029 1029			// Reserved for Automation error.
#define IGORPRO_COM_INIT_FAILED 1030			// "Initialization of COM failed."
#define IGORPRO_COM_UNEXPECTED_ERROR 1031		// "An unexpected error was occurred in an internal Automation routine."
#define IGORPRO_COM_PARAM_ERROR 1032			// "Automation client invalid parameter."
#define WAVE_USED_BY_AUTOMATION_IWave 1033						// "Wave is in use by an Automation client (IWave object)."
#define VARIABLE_USED_BY_AUTOMATION_IVariable 1034				// "Variable is in use by an Automation client (IVariable object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IEnumWaves 1035			// "Data folder is in use by an Automation client (IEnumWaves object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IWaves 1036				// "Data folder is in use by an Automation client (IWaves object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IEnumVariables 1037		// "Data folder is in use by an Automation client (IEnumVariables object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IVariables 1038			// "Data folder is in use by an Automation client (IVariables object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IDataFolder 1039			// "Data folder is in use by an Automation client (IDataFolder object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IEnumDataFolders 1040	// "Data folder is in use by an Automation client (IEnumDataFolders object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IDataFolders 1041		// "Data folder is in use by an Automation client (IDataFolders object)."
#define AUTOMATION_WAVE_KILLED_IWave 1042						// "A wave referenced by an Automation client was killed (IWave)."
#define AUTOMATION_VARIABLE_KILLED_IVariable 1043				// "A variable referenced by an Automation client was killed (IVariable)."
#define AUTOMATION_DATA_FOLDER_KILLED_IEnumWaves 1044			// "A data folder referenced by an Automation client was killed (IEnumWaves)."
#define AUTOMATION_DATA_FOLDER_KILLED_IWaves 1045				// "A data folder referenced by an Automation client was killed (IWaves)."
#define AUTOMATION_DATA_FOLDER_KILLED_IEnumVariables 1046		// "A data folder referenced by an Automation client was killed (IEnumVariables)."
#define AUTOMATION_DATA_FOLDER_KILLED_IVariables 1047			// "A data folder referenced by an Automation client was killed (IVariables)."
#define AUTOMATION_DATA_FOLDER_KILLED_IDataFolder 1048			// "A data folder referenced by an Automation client was killed (IDataFolder)."
#define AUTOMATION_DATA_FOLDER_KILLED_IEnumDataFolders 1049		// "A data folder referenced by an Automation client was killed (IEnumDataFolders)."
#define AUTOMATION_DATA_FOLDER_KILLED_IDataFolders 1050			// "A data folder referenced by an Automation client was killed (IDataFolders)."
#define AUTOMATION_RESERVED_1051 1051			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1052 1052			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1053 1053			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1054 1054			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1055 1055			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1056 1056			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1057 1057			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1058 1058			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1059 1059			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1060 1060			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1061 1061			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1062 1062			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1063 1063			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1064 1064			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1065 1065			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1066 1066			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1067 1067			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1068 1068			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1069 1069			// Reserved for Automation error.

#define SUBRANGE_DLOG_BADLABEL 1070				// "Found what appears to be a dimension label, but it is not a label in the selected wave."
#define SUBRANGE_DLOG_WRONGLABEL 1071			// "Found a dimension label, but it is the overall dimension label. You must use a label for a single element."
#define SUBRANGE_DLOG_EXPECTEDNUM 1072			// "Expected a number."
#define SUBRANGE_DLOG_MISSINGRANGE 1073			// "You must specify a range of elements for one dimension."
#define SUBRANGE_DLOG_BADEND 1074				// "Only the range dimension should specify the End."
#define SUBRANGE_DLOG_BADINCREMENT 1075			// "Only the range dimension should specify the Increment."
#define SUBRANGE_DLOG_WRONGNUMBEROFPOINTS 1076	// "Needed ^2 points in range, have ^3 points."
#define SUBRANGE_DLOG_EXPECTNUMORDIMLABEL 1077	// "Expected number or dimension label.",
#define SUBRANGE_DLOG_EXPECTNUM_GE_1 1078		// "Expected a number greater than or equal to 1.",
#define SUBRANGE_DLOG_REQUIREDPNTS_TOOBIG 1079	// "No dimension in the wave is large enough."
#define SUBRANGE_DLOG_BADRANGEDIM 1080			// "Required points too large for selected range dimension."
#define SUBRANGE_DLOG_LABEL_NOT_ALLOWED 1081	// "A dimension label is not allowed in the range dimension."
#define SUBRANGE_DLOG_INC_TOO_BIG 1082			// "This increment causes the range to excede the wave dimension size. It must be less than ^3."

#define CVODE_ILLEGAL_INPUT 1083				// "IntegrateODE reports illegal input for method 2 or 3"
#define CVODE_SETUP_FAILURE 1084				// "IntegrateODE reports setup failure for method 2 or 3"
#define CVODE_SOLVER_FAILURE 1085				// "IntegrateODE reports solver failure for method 2 or 3"
#define CVODE_BUG_ERROR 1086					// "IntegrateODE returned an unknown error code. Please contact WaveMetrics support.

#define CURVEFIT_CONF_NOALLATONCE 1087			// "It is not possible to calculate confidence bands for all-at-once fit functions."
#define DANGEROUS_FILE_COMMAND_DISABLED 1088	// "Command is disabled. It can be enabled in the Miscellaneous Setting Dialog's Misc Settings category."

#define MISMATCHED_NUMBER_OF_POINTS_ADD_BUTTON 1089		// "Different number of points for X and Y waves. Click Add button and edit range in the list below."
#define MISMATCHED_NUMBER_OF_POINTS_MORE_BUTTON 1090	// "Different number of points for X and Y waves. Click More Choices button."
#define EXPECTED_WAVE_SELECTION 1091					// "You must select a wave."

#define BAD_OBJECT_COORDINATES 1092				// Bad object coordinates.
#define CF_PLOTIT_NOEXPOFFSET 1093				// "The Graph Now button does not yet support the built-in functions exp_Xoffset or dbl_Xoffset."
#define EXPECTED_FILE_NAME 1094					// "Expected file name"	

#define CF_WRONGNUMBEROFCONSTANTS 1095			// "Wrong number of constants specified. The chosen fit function takes ^0 constants."
#define EXPECTNUMBERORAUTO 1096					// "Expected a number or \"Auto\""

#define E_FLAG_REQUIRED	1097					// "The /E flag is required for this operation."

#define INCOMPATIBLE_PACKAGE_PREFS_FILE 1098	// "The package preference file is incompatible with this version of Igor."
#define PACKAGE_PREFS_RECORD_NOT_FOUND 1099		// "Package preference record not found."
#define CORRUPT_PACKAGE_PREFS_FILE 1100			// "The package preference file is corrupt. You should delete it."
#define EXPECTED_NONNEGATIVE_INTEGER 1101		// "Expected a non-negative integer."

// End of Igor Pro 5-compatible error codes.

#define FIRST_IGOR6_ERR 1102

#define LOESS_100	1102		// "Wrong version number in lowesd.  Probably a typographic error in calling routine."
#define LOESS_101	1103		// "d>dMAX in lowesbwk (ehg131).  Need to recompile Igor with increased dimensions."
#define LOESS_102	1104		// "lowesd: liv parameter too small."
#define LOESS_103	1105		// "lowesd: lv parameter too small."
#define LOESS_104	1106		// "Fewer data values than degrees of freedom (span too small). Make /N or /SMTH value bigger."
#define LOESS_105	1107		// "k>d2MAX in l2fitHat (ehg136).  Need to recompile Igor with increased dimensions."
#define LOESS_106	1108		// "lwork too small."
#define LOESS_107	1109		// "Invalid value for kernel."
#define LOESS_108	1110		// "Invalid value for ideg."
#define LOESS_109	1111		// "lowstt only applies when kernel=1."
#define LOESS_110	1112		// "Not enough extra workspace for robustness calculation."

#define LOESS_120	1113		// "Zero-width neighborhood. Make /N or /SMTH value bigger."
#define LOESS_121	1114		// "All data on boundary of neighborhood. Make /N or /SMTH value bigger."
#define LOESS_122	1115		// "Extrapolation not allowed with blending."
#define LOESS_123	1116		// "ihat=1 (diag L) in l2fit only makes sense if z=x (eval=data)."

#define LOESS_171	1117		// "lowesd must be called first."
#define LOESS_172	1118		// "lowesf must not come between lowesb and lowese, lowesr, or lowesl."
#define LOESS_173	1119		// "lowesb must come before lowese, lowesr, or lowesl."
#define LOESS_174	1120		// "lowesb need not be called twice."
#define LOESS_175	1121		// "Need setLf=.true. for lowesl."

#define LOESS_180	1122		// "nv>nvmax in cpvert."
#define LOESS_181	1123		// "nt>20 in eval."
#define LOESS_182	1124		// "svddc failed in l2fit."
#define LOESS_183	1125		// "Didn't find edge in vleaf."
#define LOESS_184	1126		// "Zero-width cell found in vleaf."
#define LOESS_185	1127		// "Trouble descending to leaf in vleaf."
#define LOESS_186	1128		// "Insufficient workspace for lowesf."
#define LOESS_187	1129		// "Insufficient stack space."
#define LOESS_188	1130		// "lv too small for computing explicit L."

#define LOESS_191	1131		// "Computed trace L was negative; something is wrong with Loess!"
#define LOESS_192	1132		// "Computed delta was negative; something is wrong with Loess!"
#define LOESS_193	1133		// "Workspace in loread appears to be corrupted."
#define LOESS_194	1134		// "Trouble in l2fit/l2tr."
#define LOESS_195	1135		// "Only constant, linear, or quadratic local models allowed."
#define LOESS_196	1136		// "Degree must be at least 1 for vertex influence matrix."

#define LOESS_ERR1	1137		// "Specified the square of a factor predictor to be dropped when degree = 1."
#define LOESS_ERR2	1138		// "Specified the square of a predictor to be dropped with only one numeric predictor."
#define LOESS_ERR3	1139		// "Specified parametric for all predictors."

#define EXPECTED_CARET0 1140		// "Expected ^0." use BIParam1Text() to set the value of ^1.

// Stats Errors:
#define kBadMethodSpecification		1141	// "Bad method specification."
#define kInconsistentParameters		1142	// "Inconsistent choice of parameters."
#define kTooManyNaNsInData			1143	// "Too many NaNs in data."
#define kDidNotFindMissingValue		1144	// "Did not find missing value."
#define kExpectedTwoColWave			1145	// "Expected two column wave."
#define kStatsErrRes01				1146	// "AG Reserved1"
#define kStatsErrRes02				1147	// "AG Reserved2"
#define kStatsErrRes03				1148	// "AG Reserved3"
#define kStatsErrRes04				1149	// "AG Reserved4"
#define kStatsErrRes05				1150	// "AG Reserved5"
#define kStatsErrRes06				1151	// "AG Reserved6"
#define kStatsErrRes07				1152	// "AG Reserved7"
#define kStatsErrRes08				1153	// "AG Reserved8"
#define kStatsErrRes09				1154	// "AG Reserved9"
#define kStatsErrRes10				1155	// "AG Reserved10"
#define kStatsErrRes11				1156	// "AG Reserved11"
#define kStatsErrRes12				1157	// "AG Reserved12"
#define kStatsErrRes13				1158	// "AG Reserved13"
#define kStatsErrRes14				1159	// "AG Reserved14"
#define kStatsErrRes15				1160	// "AG Reserved15"
#define kStatsErrRes16				1161	// "AG Reserved16"
#define kInsufficientInput			1162	// "Insufficient input."
#define kBadInputParameter			1163	// "Inappropriate or out of range input parameter"
#define kBadWeightsData				1164	// "Bad Weights Data"
#define kPrimaryTagOverwrite		1165	// "A Tag wave can't overwrite a primary tag."

#define EXPECTED_SRCWAVE			1166	// "The srcWave parameter is required."
#define DESTFACTOR_ONE_DIM			1167	// "Destination factor waves must be one-dimensional."
#define SRCFACTOR_ONE_DIM			1168	// "Source factor waves must be one-dimensional."
#define SRCFACTORS_INAPPROPRIATE	1169	// "Source factor waves cannot be used with a two-dimensional srcWave."
#define SRCDEST_DIM_MISMATCH		1170	// "Source and destination must have the same number of dimensions."
#define LOESS_BAD_NEIGHBORS			1171	// "expected /N value to be > 1, usually about 0.1*numpnts(srcWave)."
#define LOESS_BAD_ORDER				1172	// "Expected /ORD value to be 0, 1, or 2."
#define LOESS_BAD_PASSES			1173	// "Expected /PASS value to be 1 or greater."
#define LOESS_BAD_SMOOTH			1174	// "Expected /SMTH value to be > 0.0 and <= 1.0."
#define LOESS_DestFactorNaNsNotAllowed 1175	// "NaNs not allowed in destination factor waves."
#define LOESS_SRCFACTOR_MISMATCH		1176	// "Source wave and source factor waves must have the same dimensions."
#define LOESS_DESTFACTOR_MISMATCH		1177	// "Destination wave and destination factor waves must have the same dimensions."
#define LOESS_BAD_VERBOSE				1178	// "Expected /V value to be >= 0 and <= 7."
#define LOESS_BAD_CONFIDENCEINTERVAL	1179	// "Expected /CONF value to be >= 0.0 and <= 1.0."

#define SMOOTH_BAD_DIMENSION			1180	// "Expected /DIM value to be -1, 0, 1, 2, or 3."
#define SMOOTH_BAD_THRESHOLD			1181	// "Expected /M value to be >= 0.0."
#define SMOOTH_SHORTER_THAN_WINDOW		1182	// "Length of wave ^0 can not be less than the smoothing width."
#define SMOOTH_ENDEFFECT_INCOMPATIBLE	1183	// "/E not compatible with /M."

#define EXPECTED_FILTER_COEFFICIENTS	1184	// "Expected /COEF=coefsWave, /HI, /LO, or /NMF." (FilterFIR)
#define EXPECTED_IIR_COEFFICIENTS		1185	// "Expected /COEF=coefsWave, /HI, /LO, or /N." (FilterIIR)
#define EXPECTED_IIR_COEFS_2_COLS 		1186	// "Expected a matrix wave containing 2 columns of IIR coefficient values."
#define EXPECTED_IIR_COEFS_6_COLS 		1187	// "Expected a matrix wave containing 6 columns of IIR coefficient values."

#define LFSR_MISSING_N_FLAG				1188	// "You must specify how many bits to use in the shift register (/N flag)."
#define LFSR_BADTAPNUMBER				1189	// "Bad tap specification: tap numbers must be between 1 and the number of bits in the register."
#define LFSR_ZEROINITIALCONTENTS		1190	// "The initial contents of the shift register must be non-zero (or you get nothing but zeroes forever)."

#define UNMATCHED_CONJUGATE_ZERO_OR_POLE 1191	// "Coefficients contain complex pole or zero without its matching conjugate."	
#define AMBIGUOUS_IIR_FLAGS				1192	// "/CASC/COEF/ZP is ambiguous; either supply a /COEF wave or omit one of /CASC or /ZP."

#define OLD_MAC_PRINTHANDLE_IS_OBSOLETE 1193	// "Old pre-Mac OS X page setup record is no longer supported."
#define CANT_RUN_CFM_XOP_ON_INTEL_MAC	1194	// "PowerPC XOPs can not run on Intel Macintosh when running Igor as an Intel application."
#define CANT_RUN_PPC_XOP_WITH_INTEL_IGOR	1195	// "Can't run PowerPC XOP with Intel Igor."
#define CANT_RUN_INTEL_XOP_WITH_PPC_IGOR	1196	// "Can't run Intel XOP with PowerPC Igor."
#define EXPECTED_PROC_PICT_NAME			1197		// "Expected proc pict name."
#define EXPECTED_SPECIAL_CHAR_NAME		1198		// "Expected special character name."
#define NOTEBOOK_HELPER_FILE_MISSING	1199		// "The action helper procedure file can not be found."
#define NOTEBOOK_ACTION_HELPER_BAD_PATH	1200		// "The action helper procedure file must be on the same volume as the notebook."
#define NOTEBOOK_ACTION_COMPILE_ERROR	1201		// "An error occurred during compilation of the action helper procedure file."
#define NOTEBOOK_ACTION_FILE_ALREADY_OPEN	1202	// "The action helper procedure file is already open as another type of window."
#define BAD_NOTEBOOK_ACTION_VERSION		1203	// "This action is incompatible with this version of Igor."

#define DEBUGGER_ENCOUNTERED			1204	// "Debugger statement encountered."
#define TOO_MANY_MAIN_MENUS				1205	// "There are too many main user-defined menus."
#define NOT_SVAR_OR_FCTN				1206	// "got \"^0\" instead of a string variable or string function name."

#define FUNCFIT_FITFUNC_REQUIRES_STRCFLAG 1207	// "Your fitting function takes a structure as its parameter. You must use the /STRC flag to specify an instance of the structure to be used during fitting."
#define FUNCFIT_FITFUNC_STRUCT_MISMATCH	1208	// "Your fitting funciton takes a structure of type ^0, but the /STRC flag specifies a structure of type ^1.",
#define FUNCFIT_INCOMPATIBLE_STRUCT		1209	// "The structure specified by the /STRC flag is not suitable for a fitting function."
#define FUNCFIT_FitFunctionRequestedStop 1210	// "Curve fitting stopped because the fitting function requested it."

#define EXPECT_ONE_D_WAVE				1211	// "Expected a 1-dimensional wave and got a higher dimensioned object."

#define EXPECT_CONTROL_OR_TRACE			1212	// JW 080428 for GetUserData(), "Expected name of a control or graph trace."
#define CF_XYMISMATCH					1213	// JW 080516 Curve fit was using BAD_XWave, but curve fit no longer requires matching number type. "X and Y data have different number of points."
#define JW_RESERVED_4					1214
#define JW_RESERVED_5					1215
#define JW_RESERVED_6					1216
#define JW_RESERVED_7					1217
#define CF_ONLY_ODR						1218
#define CF_ONLY_IMPLICIT				1219
#define CF_ODR_FEASIBLEONLY				1220
#define CF_ODR_BOUNDCONSTRAINTSONLY		1221
#define CF_ODR_BAD_X_RELATED_WAVE		1222
#define CF_ODR_XRESID_MISMATCH			1223
#define CF_ODR_XHOLD_MISMATCH			1224
#define CF_ODR_XWEIGHT_MISMATCH		1225
#define CF_IMPLICIT_NOMD				1226
#define CF_ODR_MISSINGCSCALWAVE			1227
#define CF_ODR_MISSINGXRESWAVE			1228
#define CF_ODR_MISSINGXHOLDWAVE			1229
#define CF_ODR_MISSINGXWGTWAVE			1230

#define LOESS_197						1231		// "NaN in predict code. Increase /N or /SMTH."
#define FCMD_AS_NOT_COMPATIBLE_WITH_SLASH_Q 1232	// "\"as\" parameter may not be used with /Q. In fact, no destination is allowed with /Q."
#define REGEX_ERROR_AT					1233		// "Regular expression error \"^1\" beginning at \"^2\"."

#define LOAD_ABORT 1234								// "File load aborted because dialogs are not allowed in threads."
#define DFREF_TYPE_INCONSISTENT	1235				// "Inconsistent type for DFREF."
#define EXPECT_DFREF_FIELD	1236					// "Expected DFREF field."
#define NOT_ALLOWED_FOR_FREE_DF 1237				//  "Can't use a free data folder in this situation."
#define NO_LOCAL_OBJ_HERE 1238						// "Can't use a free wave or an object from a free data folder in this situation."

#define NO_DFNUMERIC_WAVE_OVERWRITE 1239			// "Can't convert a Data Folder wave to or from another type."
#define NO_DF_WAVE_HERE 1240						// "Can't use a Data Folder wave in this situation."
#define NO_COMPLEX_DF_WAVES 1241					// "Can't create complex Data Folder waves."
#define NO_DF_WAVE_HERE_YET 1242					// "Can't use a Data Folder wave in this situation (yet.)"
#define EXPECTED_DF_WAVE 1243						// "Expected a Data Folder wave."
#define EXPECTED_DF_REF 1244						// "Expected a Data Folder Reference."

#define NO_WAVEREFNUMERIC_WAVE_OVERWRITE 1245		// "Can't convert a wave reference wave to or from another type."
#define NO_WAVE_WAVE_HERE 1247						// "Can't use a wave reference wave in this situation."
#define NO_WAVE_WAVE_HERE_YET 1246					// "Can't use a wave reference wave in this situation (yet.)"
#define EXPECTED_WAVE_REF 1248						// "Expected a wave reference."

#define NO_NUMERIC_WAVE_HERE	1249				// "Can't use a numeric wave in this situation."

#define WAVE_TOO_LARGE	1250						// "Wave is too large for this operation."	-- for Igor7, not used in Igor6
#define BAD_VALUE		1251						// "A provided value was not valid."
#define NUMERIC_ACCESS_ON_NON_NUMERIC_WAVE	1252	// "An attempt was made to treat a non-numeric wave (text, wave reference or DFREF) as if it were a numeric wave."
#define XOP_RETURNED_ERROR_FROM_THREAD	1253		// "An XOP returned an error from a thread."
#define LH_RES_ERR_20	1254						// "LH Reserved error 20"
#define LH_RES_ERR_21	1255						// "LH Reserved error 21"
#define LH_RES_ERR_22	1256						// "LH Reserved error 22"
#define LH_RES_ERR_23	1257						// "LH Reserved error 23"
#define LH_RES_ERR_24	1258						// "LH Reserved error 24"
#define LH_RES_ERR_25	1259						// "LH Reserved error 25"
#define LH_RES_ERR_26	1260						// "LH Reserved error 26"

#define NOT_IN_MAIN_PROC 1261						// "This construct can not be used in the main procedure window."

#define BAD_TYPE_FOR_DATE_WAVE	1262				// "Date and date/time waves must be double-precision. Use Redimension to make the wave double-precision."

#define EXPECTED_RESAMPLE_RATE	1263				// "Expected /UP, /SAME, or /RATE."
#define BAD_RATE_OR_SAME		1264				// "/RATE=rate or /SAME=wave value must be a valid number > 0."
#define BAD_OUTPUT_DIM_RATE		1265				// "DimDelta(outputWave,dim) must be a valid number > 0."

#define RESAMPLE_COEFS_WRONGPOINTS 1266				// "Needed ^2 points in /COEF=coefs, have ^3 points."

#define CFSUM_EXPECTLCBRACE_OR_STRING	1267		// "Expected left curly brace or 'string' keyword."
#define CFSUM_EXPECTSTRINGKW	1268				// "Expected 'string' keyword."
#define CFSUM_INDVAR_MISMATCH	1269				// "In your list of fit functions to be summed, at least one requires the wrong number of independent variables."

#define BAD_FREQUENCY_RATE 1270						// "Expected rate (frequency) > 0."
#define NODATA_IN_WAVE_DIM 1271						// "^1 has no data allocated in this dimension."

#define CF_NODATA		1272						// "There are no Y data points. Could be a bad sub-range, or a zero-point wave."
#define CF_NOFITCOEFS	1273						// "No fit coefficients; did you hold all of them?"
#define CF_TEXTCONSTRAINTS_NOT_THREADSAFE 1274		// "Use of /C=textWave is not allowed in a threadsafe function. Use /C={CMatrix, DVector} instead."

#define OVERWRITE_FOLDER_PERMISSION_DENIED 1275		// "You denied permission to overwrite a folder."
#define DELETE_FOLDER_PERMISSION_DENIED 1276		// "You denied permission to delete a folder."
#define FOLDER_EXISTS_NO_OVERWRITE 1277				// "The destination folder already exists and overwrite was not specified."

#define RELATIVE_INCLUDE_FROM_UNSAVED_EXP 1278		// "You must save the experiment to a file before using a relative include statement."

#define CONTEXTUAL_CAN_NOT_BE_BUILTIN_MENU 1279		// "A Menu definition with the 'contextualmenu' option cannot use the name of a built-in Igor menu."

// End of Igor Pro 6.0-compatible error codes.

#define FIRST_IGOR61_ERR 1280

#define EQUALVORONOI_REQUIRES_XYZ 1280				// "The equalVoronoiDistances= keyword applies only to XYZ contour data, yet the named contour data is a Z matrix."
#define LOESS_TIMEOUT 1281							// "Loess exceeded the allotted /TIME value. Stopping."
#define LOESS_BAD_TIME	1282						// "Expected /TIME value to be >= 0."
#define DOIGORMENU_WILLCRASH 1283					// "DoIgorMenu \"^1\", \"^2\" will crash Igor. Use Execute/P, instead."

#define CANT_OPEN_FOLDER 1284						// "Igor can not open a folder." (For when user drops an folder icon onto the Igor Pro icon in Mac OS X.)

#define CANT_INSERT_FORMATTED_TEXT_IN_PLAIN_TEXT_NOTEBOOK 1285	// "Can't insert formatted text in a plain text notebook."
#define NOTEBOOK_DATA_IS_CORRUPT_OR_OUTDATED 1286	// "The data used with the Notebook setData or Notebook zDataEnd command is corrupt or outdated."
#define BAD_NOTEBOOK_GETDATA_PARAM 1287				// "The Notebook getData parameter must be between 1 and 4."
#define CF_NO_RENTRANT_CURVEFIT 1288				// "A curve fit is being run when a curve fit is already running."

#define CANT_OVERWRITE_FOLDER_WITH_FILE 1289		// "You can't overwrite a folder with a file of the same name."
#define FILE_EXISTS_NO_OVERWRITE 1290				// "The destination file already exists and overwrite was not specified."
#define REQUIRES_LATER_OS_VERSION 1291				// "The operation could not be completed because it requires a later operating system version."

#define UNSUPPORTED_PICTURE_FORMAT 1292				// "Unknown or unsupported picture format."
#define EXPECTED_NOTEBOOK_PICTURE_NAME 1293			// "Expected name of a notebook picture."

#define EXPECTED_MODULE_NAME 1294					// "Expected module name, got \"^0\"."

#define IGOR_OBJECT_PICTURE_NOT_SELECTED 1295		// "The selection is not a single Igor-object picture."
#define IGOR_OBJECT_NOT_FOUND 1296					// "The window from which this Igor-object picture was generated no longer exists."

#define NO_SPECIAL_ACTION_WITH_DIALOG 1297			// "Can't execute a special action while a dialog is active."

#define SMOOTH_BAD_PERCENTILE 1298					// "Expected /MPCT value to be >= 0.0 and <= 100.0."

#define NB_BAD_SAVE_TYPE_FOR_SUBWINDOW 1299			// "You can save a notebook subwindow using SaveNotebook but not using /S=1 or /S=2."

#define WINDOW_MAY_NOT_BE_KILLED 1300				// "The window may not be killed at this time because a sensitive hook event is in progress."

#define BAD_CHAR_AFTER_NAME_OR_DATAFOLDER 1301		// "Bad character after name or data folder path."
#define GUIDENAMELIST_UNKNOWN_OPTION 1302			// "GuideNameList does not recognize the option \"^0\""

#define CANT_OPEN_FORMATTED_NOTEBOOK_AS_PROCWIN 1303	// "Can't open the file as a procedure file because it appears to be a formatted notebook file."

#define OPEN_MULTI_NOT_ALLOWED 1304					// "/MULT=1 is allowed only for open file dialog. Use /MULT with /R and /D=1."
#define CRVFIT_CONF_NO_ODR 1305						// "Cannot calculate confidence or prediction bands for ODR fits."
#define ODE_STOPWAVE_WRONGCOLUMNS 1306				// "The stop wave must have one column for each equation in your ODE system."
#define ODE_STOPPED_BY_STOPWAVE 1307				// "IntegrateODE stopped for a condition specified by /STOP flag."
#define ODE_PROGRAMMED_STOP 1308					// "IntegrateODE stopped by request of the derivative function."

#define RECORD_TOO_BIG_FOR_PACKED_FILE 1309			// "A record can not be written to the packed experiment file because it exceeds 4,294,967,295 bytes."
#define FLAGS_MUST_FOLLOW_FORMAT 1310				// "Flags must be placed after the format parameter" (wfprintf)

#define NAME_OR_LINE_NOT_BOTH 1311					// "ambiguous: use either /L=line OR procedure name, not both"
#define EXPECTED_NAME_OR_LINE 1312					// "need one of procedure name, /W or /L"
#define BAD_HOST 1313								// "Host window does not exist."
#define NOEXTERIOR_HOST 1314						// "Host window for an exterior panel may not be an exterior window.",

#define XOP_BAD_CALLBACK_MODE 1315					// "An XOP did a callback to Igor from a thread using an unsafe callback mode. Contact the XOP programmer."
#define XOP_CALLBACK_NOT_THREADSAFE 1316			// "An XOP called a non-threadsafe Igor callback from a thread. Contact the XOP programmer."

#define CFSUM_NULLHOLDSTRING 1317					// "In summed fit functions string, failed to get hold string."
#define CFSUM_EXPECTEDQUOTEORSVAR 1318				// "In summed fit functions hold string, expected a quoted string or name of a global string variable."

#define MAX_100_TABS 1319							// "No more than 100 tabs are allowed"

#define PERTURBATION_REQUIRES_XYZ 1320				// "The perturbation keyword applies only to XYZ contour data, yet the named contour data is a Z matrix."

#define WAVE_INDEX_OUT_OF_RANGE 1321				// "Index out of range for wave \"^0\"."
#define WAVE_NOSUCH_DIM_LABEL 1322					// "Couldn't find the given dimension item label %^1 in wave \"^0\"."

// libcurl errors.  See <http://curl.haxx.se/libcurl/c/libcurl-errors.html> for more information.
// The Igor error code constants are named the same as the CURLcode codes
// except that instead of being prefixed with CURLE_ they are prefixed with LIBCURLE_
// Many of these error codes are never actually used but they have been included
// in case they are used in the future and so that there is a complete mapping
// between the possible libcurl error codes and the Igor error codes.
#define LIBCURLE_UNKNOWN 1323						/* "The network library returned an unknown error." */
#define LIBCURLE_UNSUPPORTED_PROTOCOL 1324			/* "The specified URL uses an unsupported protocol." */
#define LIBCURLE_FAILED_INIT 1325					/* "The network library failed to properly initialize." */
#define LIBCURLE_URL_MALFORMAT 1326					/* "The specified URL is missing or not properly formatted." */
#define LIBCURLE_COULDNT_RESOLVE_PROXY 1327			/* "The specified proxy host could not be resolved." */
#define LIBCURLE_COULDNT_RESOLVE_HOST 1328			/* "The specified host could not be resolved." */
#define LIBCURLE_COULDNT_CONNECT 1329				/* "Could not connect to server." */
#define LIBCURLE_FTP_WEIRD_SERVER_REPLY 1330		/* "FTP: The server gave an unexpected reply." */
#define LIBCURLE_REMOTE_ACCESS_DENIED 1331			/* "Access to remote resource was denied." */
#define LIBCURLE_FTP_WEIRD_PASS_REPLY 1332			/* "FTP: The server gave an unexpected reply to the password command." */
#define LIBCURLE_FTP_WEIRD_PASV_REPLY 1333			/* "FTP: The server gave an unexpected reply to the PASV command." */
#define LIBCURLE_FTP_WEIRD_227_FORMAT 1334			/* "FTP: Got an unexpected response format to the PASV command." */
#define LIBCURLE_FTP_CANT_GET_HOST 1335				/* "FTP: The host supplied by the PASV response could not be looked up." */
#define LIBCURLE_FTP_COULDNT_SET_TYPE 1336			/* "FTP: Could not set the file type." */
#define LIBCURLE_PARTIAL_FILE 1337					/* "A file transfer was shorter or larger than expected." */
#define LIBCURLE_FTP_COULDNT_RETR_FILE 1338			/* "FTP: Could not retrieve the specified file." */
#define LIBCURLE_QUOTE_ERROR 1339					/* "A custom command returned an error code." */
#define LIBCURLE_HTTP_RETURNED_ERROR 1340			/* "HTTP: The server returned an error code." */
#define LIBCURLE_WRITE_ERROR 1341					/* "An error occurred when writing received data to a local file." */
#define LIBCURLE_UPLOAD_FAILED 1342					/* "Failure in starting the upload." */
#define LIBCURLE_READ_ERROR 1343					/* "There was a problem reading a local file." */
#define LIBCURLE_OUT_OF_MEMORY 1344					/* "The network library is out of memory." */
#define LIBCURLE_OPERATION_TIMEDOUT 1345			/* "The network library timeout was reached." */
#define LIBCURLE_FTP_PORT_FAILED 1346				/* "FTP: The PORT command failed." */
#define LIBCURLE_FTP_COULDNT_USE_REST 1347			/* "FTP: The REST command failed." */
#define LIBCURLE_RANGE_ERROR 1348					/* "The server does not support or accept range requests. " */
#define LIBCURLE_HTTP_POST_ERROR 1349				/* "HTTP: There was an error in executing the POST method." */
#define LIBCURLE_SSL_CONNECT_ERROR 1350				/* "A problem occurred somewhere in the SSL/TLS handshake." */
#define LIBCURLE_BAD_DOWNLOAD_RESUME 1351			/* "The download could not be resumed because the specified offset was out of the file boundary." */
#define LIBCURLE_FILE_COULDNT_READ_FILE 1352		/* "Could not open the specified file using the file:// scheme." */
#define LIBCURLE_LDAP_CANNOT_BIND 1353				/* "LDAP: Bind operation failed." */
#define LIBCURLE_LDAP_SEARCH_FAILED 1354			/* "LDAP: Search failed." */
#define LIBCURLE_FUNCTION_NOT_FOUND 1355			/* "A required function in the library was not found." */
#define LIBCURLE_ABORTED_BY_CALLBACK 1356			/* "Operation was aborted by an application callback." */
#define LIBCURLE_BAD_FUNCTION_ARGUMENT 1357			/* "A network function was called with a bad parameter." */
#define LIBCURLE_INTERFACE_FAILED 1358				/* "A specified outgoing interface could not be used." */
#define LIBCURLE_TOO_MANY_REDIRECTS 1359			/* "Number of redirects hit maximum amount." */
#define LIBCURLE_UNKNOWN_TELNET_OPTION 1360			/* "An unknown telnet option was specified." */
#define LIBCURLE_TELNET_OPTION_SYNTAX 1361			/* "A telnet option string was illegally formatted." */
#define LIBCURLE_PEER_FAILED_VERIFICATION 1362		/* "SSL peer certificate or SSH remote key was not OK." */
#define LIBCURLE_GOT_NOTHING 1363					/* "The server response was empty." */
#define LIBCURLE_SSL_ENGINE_NOTFOUND 1364			/* "The specified SSL cryptographic engine was not found." */
#define LIBCURLE_SSL_ENGINE_SETFAILED 1365			/* "Can not set the specified SSL crypto engine as default." */
#define LIBCURLE_SEND_ERROR 1366					/* "Failed sending network data." */
#define LIBCURLE_RECV_ERROR 1367					/* "Failed receiving network data." */
#define LIBCURLE_SSL_CERTPROBLEM 1368				/* "Problem with the local SSL certificate." */
#define LIBCURLE_SSL_CIPHER 1369					/* "Couldn't use specified SSL cipher." */
#define LIBCURLE_SSL_CACERT 1370					/* "Peer certificate cannot be authenticated with known CA certificates." */
#define LIBCURLE_BAD_CONTENT_ENCODING 1371			/* "Unrecognized transfer encoding." */
#define LIBCURLE_LDAP_INVALID_URL 1372				/* "LDAP: Invalid URL." */
#define LIBCURLE_FILESIZE_EXCEEDED 1373				/* "Maximum file size exceeded." */
#define LIBCURLE_USE_SSL_FAILED 1374				/* "FTP: Requested SSL level failed." */
#define LIBCURLE_SEND_FAIL_REWIND 1375				/* "Send failed since rewinding of the data stream failed." */
#define LIBCURLE_SSL_ENGINE_INITFAILED 1376			/* "Failure to initialize the specified SSL cryptographic engine." */
#define LIBCURLE_LOGIN_DENIED 1377					/* "Login failure." */
#define LIBCURLE_TFTP_NOTFOUND 1378					/* "TFTP: File not found on server." */
#define LIBCURLE_TFTP_PERM 1379						/* "TFTP: Permission problem on server." */
#define LIBCURLE_REMOTE_DISK_FULL 1380				/* "Out of disk space on the remote server." */
#define LIBCURLE_TFTP_ILLEGAL 1381					/* "TFTP: Illegal operation." */
#define LIBCURLE_TFTP_UNKNOWNID 1382				/* "TFTP: Unknown transfer ID." */
#define LIBCURLE_REMOTE_FILE_EXISTS 1383			/* "Remote file already exists and will not be overwritten." */
#define LIBCURLE_TFTP_NOSUCHUSER 1384				/* "TFTP:  No such user exists." */
#define LIBCURLE_CONV_FAILED 1385					/* "Character conversion failed." */
#define LIBCURLE_CONV_REQD 1386						/* "Caller must register conversion callbacks." */
#define LIBCURLE_SSL_CACERT_BADFILE 1387			/* "Problem reading the SSL CA certificate." */
#define LIBCURLE_REMOTE_FILE_NOT_FOUND 1388			/* "The resource referenced in the URL does not exist." */
#define LIBCURLE_SSH 1389							/* "SSH: An unspecified error occurred during the session." */
#define LIBCURLE_SSL_SHUTDOWN_FAILED 1390			/* "Failed to shut down the SSL connection." */
#define LIBCURLE_AGAIN 1391							/* "Socket not ready for send/recv." */
#define LIBCURLE_SSL_CRL_BADFILE 1392				/* "Failed to load CRL file." */
#define LIBCURLE_SSL_ISSUER_ERROR 1393				/* "Issuer check against peer certificate failed." */
#define LIBCURLE_FTP_PRET_FAILED 1394				/* "FTP: The server did not accept the PRET command." */
#define LIBCURLE_RTSP_CSEQ_ERROR 1395				/* "RTSP: CSeq mismatch or invalid CSeq." */
#define LIBCURLE_RTSP_SESSION_ERROR 1396			/* "RTSP: Mismatch of session identifiers." */
#define LIBCURLE_FTP_BAD_FILE_LIST 1397				/* "FTP: Unable to parse file list." */
#define LIBCURLE_CHUNK_FAILED 1398					/* "Chunk callback reported error." */
#define LIBCURLE_RESERVED_1399 1399					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1400 1400					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1401 1401					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1402 1402					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1403 1403					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1404 1404					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1405 1405					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1406 1406					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1407 1407					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1408 1408					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1409 1409					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1410 1410					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1411 1411					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1412 1412					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1413 1413					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1414 1414					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1415 1415					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1416 1416					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1417 1417					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1418 1418					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1419 1419					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1420 1420					/* "Reserved for future libcurl error." */
#define LIBCURLE_RESERVED_1421 1421					/* "Reserved for future libcurl error." */
// end of libcurl errors

#define NULL_DATAFOLDER_OP 1422						// "Attempt to operate on a NULL data folder reference."
#define OH_EXPECTED_DATAFOLDER_REF 1423				// "Expected a data folder reference."

// From WM.h - More error codes that Igor understands
#define FRST_WM_ERR 1500				/* ID of first WM error string */
#define LAST_WM_ERR 1999				/* ID of last POSSIBLE WM error string */
enum {
	CANT_OVERWRITE_OPEN_FILE = FRST_WM_ERR,	/* "can't overwrite open file" */
	CANT_HANDLE_FILE_TYPE,				/* "can't handle this type of file" */
	FILE_VERSION_TOO_OLD,				/* "this file version is no longer supported" */
	FILE_VERSION_TOO_NEW,				/* "this file version requires a later version of the program" */
	FILE_IS_READ_ONLY,					/* "the file is open for reading only" */
	FILE_IS_OPEN_SOMEWHERE_ELSE,		/* "the file Ò^0Ó could not be written because there is another copy of the file open" */
	GENERAL_BAD_VIBS,					// "Possible bug (or corrupt data, or ...)"
	EM_SECTION_NOT_FOUND,				// "couldn't locate a publisher(subscriber)"
	EM_NO_SUCH_WIN,						// "couldn't find a window needed by a publisher(subscriber)"
	CORRUPT_EDITION_INFO,				// "file is damaged -- corrupt publisher(subscriber) info"
	EM_OLD_VERSION,						// "obsolete version"
	DETACH_RESOURCE_FAILED,				// "an attempt to detach a resource failed"
	WM_BAD_FILE_NAME,					// "ill-formed file name"
	WM_FILE_NAME_TOO_LONG,				// "a file name is limited to 255 characters"
	DIRECTORY_REFERENCED_TWICE,			// "an alias has created two references to the same folder"
	BALLOON_INSERT_ERR,					// "an error occurred while creating a help balloon"
	MINUS_ONE_ERR,						// "an error occurred"
	WM_FILEPATH_TOO_LONG,				// "the file and path is too long"
	CANT_WRITE_LOCKED_FILE,				// "the file Ò^0Ó cannot be overwritten because it is locked"
	WM_PATH_TOO_LONG,					// "the path to a file or folder is too long"
	WM_BAD_VOLUME_SPEC,					// "the specified volume can not be found"
	WM_BAD_DIR_SPEC,					// "the specified directory can not be found"
	WM_BAD_FILE_OR_FOLDER_SPEC,			// "the file or folder can not be found"
	WM_BAD_DIRID,						// "BUG: A bad directory ID was used"
	WM_BAD_WDREFNUM,					// "BUG: A bad working directory reference number was used"
	WM_BAD_FILE_REFNUM,					// "BUG: A bad file reference number was used"
	WM_DIRNAME_TOO_LONG,				// "a directory name is limited to 255 characters"
	WM_FILENAME_TOO_LONG,				// "a file name is limited to 255 characters"
	WM_FILE_TOO_LARGE,					// "Can't handle a file with more than 2^32 bytes"
	WM_UNKNOWN_ERROR,					// "An error of an unknown nature occurred"
	WM_PRINTER_DRIVER_NOT_OPEN,			// "The printer driver is not open"
	WM_BAD_PATH_SYNTAX,					// "The path is not properly formed"
	WM_DATEFORMAT_BAD_SYSTEMDATE,		// "Expected a system date format code between 0 and 2."
	WM_DATEFORMAT_BAD_LANGUAGE,			// "Expected a language code between 1 and 16"
	WM_DATEFORMAT_BAD_YEARFORMAT,		// "Expected a year format code betwen 1 and 2."
	WM_DATEFORMAT_BAD_MONTHFORMAT,		// "Expected a month format between 1 and 2".
	WM_DATEFORMAT_BAD_DAYOFMONTHFORMAT,	// "Expected a day-of-month format between 1 and 2."
	WM_DATEFORMAT_BAD_DAYOFWEEKFORMAT,	// "Expected a day-of-month format between 1 and 2."
	WM_DATEFORMAT_BAD_LAYOUT,			// "The layout of the date format is incorrect."
	WM_DATEFORMAT_BAD_SEPARATOR,		// "A date separator is limited to 15 characters."
	WM_DATEFORMAT_BAD_PIVOTYEAR,		// "Expected a pivot year between 4 and 40."
	WM_DATEFORMAT_BAD_VERSION,			// "Unrecognized date format."
	EXPECTED_WM_LANGUAGE_NAME,			// "Expected the name of a supported language."
	EXPECTED_DATE,						// "Expected date."
	DUPLICATE_EXPECTED_TIME,			// "Expected time."	// This is a duplicate of EXPECTED_TIME. Use EXPECTED_TIME instead.
	EXPECTED_DATETIME,					// "Expected date/time (date<space>time)."
	URL_MISSING_BRACKETS,				// "The URL must be enclosed in angle brackets (e.g., <http://www.wavemetrics.com>)."
	URL_TRAILING_JUNK,					// "There are extraneous characters in the URL after the trailing angle bracket."
	URL_MISSING_PROTOCOL,				// "The URL must include \"http://\" or \"ftp://\"."
	URL_MISSING_HTTP,					// "The URL must include \"http://\" (e.g., <http://www.wavemetrics.com>)."
	URL_MISSING_FTP,					// "The URL must include \"ftp://\" (e.g., <ftp://ftp.wavemetrics.com>)."
	URL_TOO_LONG,						// "The URL exceeds 255 characters in length."
	WM_FTP_NOT_AVAILABLE,				// "FTP support is not available."	// Should never happen.
	WM_FTP_HOST_ERROR,					// "The FTP host returned an unknown error."
	WM_FTP_OPERATION_TIMED_OUT,			// "The operation has timed out."
	WM_FTP_SERVER_ERROR,				// "The FTP server returned an error."
	WM_INVALID_URL,						// "Invalid URL."
	WM_FTP_INCORRECT_USERNAME,			// "Incorrect user name."
	WM_FTP_INCORRECT_PASSWORD,			// "Incorrect password."
	WM_FTP_CANT_LOG_IN,					// "The request to connect to and log in to FTP server failed."
	WM_FTP_ITEM_NOT_FOUND,				// "The item requested via FTP could not be found."
	WM_FTP_CANT_CONNECT,				// "The attempt to connect to the FTP server failed."
	WM_FTP_BAD_DIAGNOSTIC_MODE,			// "Expected an FTP diagnostic mode code from 0 to 7."
	WM_FTP_BAD_TRANSFER_TYPE,			// "Expected an FTP transfer type code from 0 (binary) to 1 (ASCII)."
	WM_FTP_BAD_ENCODING,				// "Expected an FTP encoding code from 0 (no binhex) to 1 (binhex)."
	WM_FTP_BAD_PROGRESS,				// "Expected an FTP progress code from 0 (off) to 1 (on)."
	EXPECTED_STRING,					// "Expected string"
	EXPECTED_OCTAL_NUMBER,				// "Expected octal number"
	EXPECTED_HEX_NUMBER,				// "Expected hex number"
	EXPECTED_TIME_OF_DAY,				// "Expected time of day in hh:mm:ss format"
	WM_BAD_NUMTYPE,						// "Invalid or unsupported number type"
	PRINTER_NOT_FOUND,					// "Printer not found"
	UNKNOWN_SPECIAL_DIRECTORY,			// "Unknown special directory"
	SAFE_SAVE_CANT_EXCHANGE_FILES,		// "Safe save can't finalize the save"
	SAFE_SAVE_CANT_DELETE_TEMP_FILE,	// "Safe save can't delete temporary file"
	SAFE_SAVE_CANT_DELETE_ORIGINAL_FILE,// "Safe save can't delete original file"
	SAFE_SAVE_CANT_RENAME_NEW_FILE,		// "Safe save can't rename new file"
	SAFE_SAVE_CANT_RENAME_TEMP_FILE,	// "Safe save can't rename temporary file"
	CANT_CONVERT_FILE_PATH_TO_UTF2,		// "Can't convert a system file path to a Unicode file path"
	CANT_CONVERT_UTF2_TO_FILE_PATH,		// "Can't convert a Unicode file path to a system file path"
};

#endif // IGOR_ERRORS_H
