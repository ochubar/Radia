
#ifndef MOSMESA_H
#define MOSMESA_H

#include <GL/gl.h>

#ifndef GLTHREAD_H
#include "glthread.h"
#endif

//-------------------------------------------------------------------------

/* Create a macro so that asm functions can be linked into compilers other
 * than GNU C
 */
#ifndef _ASMAPI
#if defined(WIN32) && !defined(BUILD_FOR_SNAP)/* was: !defined( __GNUC__ ) && !defined( VMS ) && !defined( __INTEL_COMPILER )*/
#define _ASMAPI __cdecl
#else
#define _ASMAPI
#endif
#ifdef	PTR_DECL_IN_FRONT
#define	_ASMAPIP * _ASMAPI
#else
#define	_ASMAPIP _ASMAPI *
#endif
#endif

//-------------------------------------------------------------------------
/*
 * Forward declaration of display list data types:
 */
union Mnode;
typedef union Mnode MNode;

#define FEATURE_ARB_vertex_buffer_object 1

#define OSMESA_RGBA GL_RGBA
#define MMAX_TEXTURE_COORD_UNITS 8
#define MMAX_PROGRAM_MATRICES 8

typedef struct Mosmesa_context *MOSMesaContext;

typedef struct M__GLcontextRec MGLcontext;
typedef struct M__GLcontextModesRec MGLvisual;
typedef struct Mgl_frame_buffer MGLframebuffer;

#define MMAT_ATTRIB_MAX                     12
/** Maximum attribute stack depth */
#define MMAX_ATTRIB_STACK_DEPTH 16
/** Number of texture units - GL_ARB_multitexture */
#define MMAX_TEXTURE_UNITS 8
/** Maximum number of lights */
#define MMAX_LIGHTS 8
/** Maximum pixel map lookup table size */
#define MMAX_PIXEL_MAP_TABLE 256
/** Number of texture units - GL_ARB_multitexture */
#define MMAX_TEXTURE_UNITS 8
/** Number of 1D/2D texture mipmap levels */
#define MMAX_TEXTURE_LEVELS 12
/** Maximum user-defined clipping planes */
#define MMAX_CLIP_PLANES 6
/** Size of histogram tables */
#define MHISTOGRAM_TABLE_SIZE 256
/** Max convolution filter width */
#define MMAX_CONVOLUTION_WIDTH 9
/** Max convolution filter height */
#define MMAX_CONVOLUTION_HEIGHT 9
/** Maximum client attribute stack depth */
#define MMAX_CLIENT_ATTRIB_STACK_DEPTH 16
/**
 * \name Separate numbers of texture coordinates and texture image units.
 *
 * These values will eventually replace most instances of MAX_TEXTURE_UNITS.
 * We should always have MAX_TEXTURE_COORD_UNITS <= MAX_TEXTURE_IMAGE_UNITS.
 * And, GL_MAX_TEXTURE_UNITS <= MAX_TEXTURE_COORD_UNITS.
 */
/*@{*/
#define MMAX_TEXTURE_COORD_UNITS 8
#define MMAX_TEXTURE_IMAGE_UNITS 8
/*@}*/
/** Maximum Name stack depth */
#define MMAX_NAME_STACK_DEPTH 64
/** Maximum Number of auxillary color buffers */
#define MMAX_AUX_BUFFERS 4

#define MMAX_NV_FRAGMENT_PROGRAM_TEMPS         96
#define MMAX_NV_FRAGMENT_PROGRAM_INPUTS        12
#define MMAX_NV_FRAGMENT_PROGRAM_OUTPUTS        3
#define MMAX_NV_FRAGMENT_PROGRAM_PARAMS        64

#define MMAX_NV_VERTEX_PROGRAM_TEMPS         12
#define MMAX_NV_VERTEX_PROGRAM_PARAMS        96
#define MMAX_NV_VERTEX_PROGRAM_INPUTS        16
#define MMAX_NV_VERTEX_PROGRAM_OUTPUTS       15

#define MMAX_PROGRAM_LOCAL_PARAMS 96
#define MMAX_TEXTURE_IMAGE_UNITS 8

#define MTABLE_SIZE 1023  /**< Size of lookup table/array */
#define MSHINE_TABLE_SIZE 256	/**< Material shininess lookup table sizes */
#define MMAX_DLIST_EXT_OPCODES 16

#define MMAX_PIPELINE_STAGES     30

#define MEXP_TABLE_SIZE 512	/**< Specular exponent lookup table sizes */

#define MMAX_FACES  6

#define MACCUM_BITS 8
#define MSTENCIL_BITS 8
#define MCHAN_BITS 8

/** 
 * Maximum viewport/image width. Must accomodate all texture sizes too. 
 */
#define MMAX_WIDTH 4096
/** Maximum viewport/image height */
#define MMAX_HEIGHT 4096

#define M_TNL_MAX_ATTR_CODEGEN 16 
#define MVERT_BUFFER_SIZE 2048	/* 8kbytes */
#define MTNL_MAX_PRIM 16
#define MTNL_MAX_COPIED_VERTS 3

/**
 * Accumulation buffer data type.
 */
#if MACCUM_BITS==8
   typedef GLbyte MGLaccum;
#elif MACCUM_BITS==16
   typedef GLshort MGLaccum;
#elif MACCUM_BITS==32
   typedef GLfloat MGLaccum;
#else
//#  error "illegal number of accumulation bits"
#endif

/**
 * Stencil buffer data type.
 */
#if MSTENCIL_BITS==8
   typedef GLubyte MGLstencil;
#  define MSTENCIL_MAX 0xff
#elif MSTENCIL_BITS==16
   typedef GLushort MGLstencil;
#  define MSTENCIL_MAX 0xffff
#else
//#error "illegal number of stencil bits"
#endif

/**
 * Color channel data type.
 */
#if MCHAN_BITS == 8
   typedef GLubyte MGLchan;
#define MCHAN_MAX 255
#define MCHAN_MAXF 255.0F
#define MCHAN_TYPE GL_UNSIGNED_BYTE
#elif MCHAN_BITS == 16
   typedef GLushort MGLchan;
#define MCHAN_MAX 65535
#define MCHAN_MAXF 65535.0F
#define MCHAN_TYPE GL_UNSIGNED_SHORT
#elif MCHAN_BITS == 32
   typedef GLfloat MGLchan;
#define MCHAN_MAX 1.0
#define MCHAN_MAXF 1.0F
#define MCHAN_TYPE GL_FLOAT
#else
#error "illegal number of color channel bits"
#endif

/**
 * Depth buffer data type.
 *
 * \note Must be 32-bits!
 */
typedef GLuint MGLdepth;  

/**
 * Fixed point data type.
 */
typedef int MGLfixed;

/**
 * These define the aliases between numbered vertex attributes and
 * conventional OpenGL vertex attributes.  We use these values in
 * quite a few places.  
 *
 * New in Mesa 4.1.
 */
enum {
	MVERT_ATTRIB_POS = 0,
	MVERT_ATTRIB_WEIGHT = 1,
	MVERT_ATTRIB_NORMAL = 2,
	MVERT_ATTRIB_COLOR0 = 3,
	MVERT_ATTRIB_COLOR1 = 4,
	MVERT_ATTRIB_FOG = 5,
	MVERT_ATTRIB_SIX = 6,
	MVERT_ATTRIB_SEVEN = 7,
	MVERT_ATTRIB_TEX0 = 8,
	MVERT_ATTRIB_TEX1 = 9,
	MVERT_ATTRIB_TEX2 = 10,
	MVERT_ATTRIB_TEX3 = 11,
	MVERT_ATTRIB_TEX4 = 12,
	MVERT_ATTRIB_TEX5 = 13,
	MVERT_ATTRIB_TEX6 = 14,
	MVERT_ATTRIB_TEX7 = 15,
	MVERT_ATTRIB_MAX = 16
};

/**
 * Accumulation buffer attributes.
 */
struct Mgl_accum_attrib {
   GLfloat ClearColor[4];	/**< Accumulation buffer clear color */
};

/**
 * Color buffers attributes.
 */
struct Mgl_colorbuffer_attrib {
   GLuint ClearIndex;			/**< Index to use for glClear */
   GLclampf ClearColor[4];		/**< Color to use for glClear */

   GLuint IndexMask;			/**< Color index write mask */
   GLubyte ColorMask[4];		/**< Each flag is 0xff or 0x0 */

   GLenum DrawBuffer;			/**< Which buffer to draw into */
   GLbitfield _DrawDestMask;		/**< bitmask of DD_*_BIT bits */

   /** 
    * \name alpha testing
    */
   /*@{*/
   GLboolean AlphaEnabled;		/**< Alpha test enabled flag */
   GLenum AlphaFunc;			/**< Alpha test function */
   GLclampf AlphaRef;			/**< Alpha reference value */
   /*@}*/

   /** 
    * \name Blending
    */
   /*@{*/
   GLboolean BlendEnabled;		/**< Blending enabled flag */
   GLenum BlendSrcRGB;			/**< Blending source operator */
   GLenum BlendDstRGB;			/**< Blending destination operator */
   GLenum BlendSrcA;			/**< GL_INGR_blend_func_separate */
   GLenum BlendDstA;			/**< GL_INGR_blend_func_separate */
   GLenum BlendEquationRGB;		/**< Blending equation */
   GLenum BlendEquationA;		/**< GL_EXT_blend_equation_separate */
   GLfloat BlendColor[4];		/**< Blending color */
   /*@}*/

   /** 
    * \name Logic op
    */
   /*@{*/
   GLenum LogicOp;			/**< Logic operator */
   GLboolean IndexLogicOpEnabled;	/**< Color index logic op enabled flag */
   GLboolean ColorLogicOpEnabled;	/**< RGBA logic op enabled flag */
   GLboolean _LogicOpEnabled;		/**< RGBA logic op + EXT_blend_logic_op enabled flag */
   /*@}*/

   GLboolean DitherFlag;		/**< Dither enable flag */
};

/**
 * Current attributes.
 */
struct Mgl_current_attrib {
   /**
    * \name Values valid only when FLUSH_VERTICES has been called.
    */
   /*@{*/
   GLfloat Attrib[MVERT_ATTRIB_MAX][4];		/**< Current vertex attributes
						  *  indexed by VERT_ATTRIB_* */
   GLfloat Index;				/**< Current color index */
   GLboolean EdgeFlag;				/**< Current edge flag */
   /*@}*/

   /**
    * \name Values are always valid.  
    * 
    * \note BTW, note how similar this set of attributes is to the SWvertex
    * data type in the software rasterizer...
    */
   /*@{*/
   GLfloat RasterPos[4];			/**< Current raster position */
   GLfloat RasterDistance;			/**< Current raster distance */
   GLfloat RasterColor[4];			/**< Current raster color */
   GLfloat RasterSecondaryColor[4];             /**< Current raster secondary color */
   GLfloat RasterIndex;				/**< Current raster index */
   GLfloat RasterTexCoords[MMAX_TEXTURE_UNITS][4];/**< Current raster texcoords */
   GLboolean RasterPosValid;			/**< Raster pos valid flag */
   /*@}*/
};

/**
 * Depth buffer attributes.
 */
struct Mgl_depthbuffer_attrib {
   GLenum Func;			/**< Function for depth buffer compare */
   GLclampd Clear;		/**< Value to clear depth buffer to */
   GLboolean Test;		/**< Depth buffering enabled flag */
   GLboolean Mask;		/**< Depth buffer writable? */
   GLboolean OcclusionTest;	/**< GL_HP_occlusion_test */
   GLboolean BoundsTest;        /**< GL_EXT_depth_bounds_test */
   GLfloat BoundsMin, BoundsMax;/**< GL_EXT_depth_bounds_test */
};

/** 
 * Hint attributes.
 * 
 * Values are always one of GL_FASTEST, GL_NICEST, or GL_DONT_CARE.
 */
struct Mgl_hint_attrib {
   GLenum PerspectiveCorrection;
   GLenum PointSmooth;
   GLenum LineSmooth;
   GLenum PolygonSmooth;
   GLenum Fog;
   GLenum ClipVolumeClipping;   /**< GL_EXT_clip_volume_hint */
   GLenum TextureCompression;   /**< GL_ARB_texture_compression */
   GLenum GenerateMipmap;       /**< GL_SGIS_generate_mipmap */
};

/**
 * Eval attributes.
 */
struct Mgl_eval_attrib {
   /**
    * \name Enable bits 
    */
   /*@{*/
   GLboolean Map1Color4;
   GLboolean Map1Index;
   GLboolean Map1Normal;
   GLboolean Map1TextureCoord1;
   GLboolean Map1TextureCoord2;
   GLboolean Map1TextureCoord3;
   GLboolean Map1TextureCoord4;
   GLboolean Map1Vertex3;
   GLboolean Map1Vertex4;
   GLboolean Map1Attrib[16];  /* GL_NV_vertex_program */
   GLboolean Map2Color4;
   GLboolean Map2Index;
   GLboolean Map2Normal;
   GLboolean Map2TextureCoord1;
   GLboolean Map2TextureCoord2;
   GLboolean Map2TextureCoord3;
   GLboolean Map2TextureCoord4;
   GLboolean Map2Vertex3;
   GLboolean Map2Vertex4;
   GLboolean Map2Attrib[16];  /* GL_NV_vertex_program */
   GLboolean AutoNormal;
   /*@}*/
   
   /**
    * \name Map Grid endpoints and divisions and calculated du values
    */
   /*@{*/
   GLint MapGrid1un;
   GLfloat MapGrid1u1, MapGrid1u2, MapGrid1du;
   GLint MapGrid2un, MapGrid2vn;
   GLfloat MapGrid2u1, MapGrid2u2, MapGrid2du;
   GLfloat MapGrid2v1, MapGrid2v2, MapGrid2dv;
   /*@}*/
};

/**
 * Fog attributes.
 */
struct Mgl_fog_attrib {
   GLboolean Enabled;		/**< Fog enabled flag */
   GLfloat Color[4];		/**< Fog color */
   GLfloat Density;		/**< Density >= 0.0 */
   GLfloat Start;		/**< Start distance in eye coords */
   GLfloat End;			/**< End distance in eye coords */
   GLfloat Index;		/**< Fog index */
   GLenum Mode;			/**< Fog mode */
   GLboolean ColorSumEnabled;
   GLenum FogCoordinateSource;  /**< GL_EXT_fog_coord */
};

/**
 * Transformation attributes.
 */
struct Mgl_transform_attrib {
   GLenum MatrixMode;				/**< Matrix mode */
   GLfloat EyeUserPlane[MMAX_CLIP_PLANES][4];	/**< User clip planes */
   GLfloat _ClipUserPlane[MMAX_CLIP_PLANES][4];	/**< derived */
   GLuint ClipPlanesEnabled;                    /**< on/off bitmask */
   GLboolean Normalize;				/**< Normalize all normals? */
   GLboolean RescaleNormals;			/**< GL_EXT_rescale_normal */
   GLboolean RasterPositionUnclipped;           /**< GL_IBM_rasterpos_clip */

   GLboolean CullVertexFlag;	/**< True if GL_CULL_VERTEX_EXT is enabled */
   GLfloat CullEyePos[4];
   GLfloat CullObjPos[4];
};

/**
 * Light.
 */
struct Mgl_light {
   struct Mgl_light *next;	/**< double linked list with sentinel */
   struct Mgl_light *prev;

   GLfloat Ambient[4];		/**< ambient color */
   GLfloat Diffuse[4];		/**< diffuse color */
   GLfloat Specular[4];		/**< specular color */
   GLfloat EyePosition[4];	/**< position in eye coordinates */
   GLfloat EyeDirection[4];	/**< spotlight dir in eye coordinates */
   GLfloat SpotExponent;
   GLfloat SpotCutoff;		/**< in degrees */
   GLfloat _CosCutoff;		/**< = MAX(0, cos(SpotCutoff)) */
   GLfloat ConstantAttenuation;
   GLfloat LinearAttenuation;
   GLfloat QuadraticAttenuation;
   GLboolean Enabled;		/**< On/off flag */

   /** 
    * \name Derived fields
    */
   /*@{*/
   GLuint _Flags;		/**< State */

   GLfloat _Position[4];	/**< position in eye/obj coordinates */
   GLfloat _VP_inf_norm[3];	/**< Norm direction to infinite light */
   GLfloat _h_inf_norm[3];	/**< Norm( _VP_inf_norm + <0,0,1> ) */
   GLfloat _NormDirection[4];	/**< normalized spotlight direction */
   GLfloat _VP_inf_spot_attenuation;

   GLfloat _SpotExpTable[MEXP_TABLE_SIZE][2];  /**< to replace a pow() call */
   GLfloat _MatAmbient[2][3];	/**< material ambient * light ambient */
   GLfloat _MatDiffuse[2][3];	/**< material diffuse * light diffuse */
   GLfloat _MatSpecular[2][3];	/**< material spec * light specular */
   GLfloat _dli;		/**< CI diffuse light intensity */
   GLfloat _sli;		/**< CI specular light intensity */
   /*@}*/
};

/**
 * Light model.
 */
struct Mgl_lightmodel {
   GLfloat Ambient[4];		/**< ambient color */
   GLboolean LocalViewer;	/**< Local (or infinite) view point? */
   GLboolean TwoSide;		/**< Two (or one) sided lighting? */
   GLenum ColorControl;		/**< either GL_SINGLE_COLOR
				 *    or GL_SEPARATE_SPECULAR_COLOR */
};

/**
 * Material.
 */
struct Mgl_material
{
   GLfloat Attrib[MMAT_ATTRIB_MAX][4];
};

/**
 * Lighting attributes.
 */
struct Mgl_light_attrib {
   struct Mgl_light Light[MMAX_LIGHTS];	/**< Array of lights */
   struct Mgl_lightmodel Model;		/**< Lighting model */

   /**
    * Must flush FLUSH_VERTICES before referencing:
    */
   /*@{*/
   struct Mgl_material Material; 	/**< Includes front & back values */
   /*@}*/

   GLboolean Enabled;			/**< Lighting enabled flag */
   GLenum ShadeModel;			/**< GL_FLAT or GL_SMOOTH */
   GLenum ColorMaterialFace;		/**< GL_FRONT, BACK or FRONT_AND_BACK */
   GLenum ColorMaterialMode;		/**< GL_AMBIENT, GL_DIFFUSE, etc */
   GLuint ColorMaterialBitmask;		/**< bitmask formed from Face and Mode */
   GLboolean ColorMaterialEnabled;

   struct Mgl_light EnabledList;         /**< List sentinel */

   /** 
    * Derived for optimizations: 
    */
   /*@{*/
   GLboolean _NeedEyeCoords;		
   GLboolean _NeedVertices;		/**< Use fast shader? */
   GLuint  _Flags;		        /**< LIGHT_* flags, see above */
   GLfloat _BaseColor[2][3];
   /*@}*/
};

/**
 * Line attributes.
 */
struct Mgl_line_attrib {
   GLboolean SmoothFlag;	/**< GL_LINE_SMOOTH enabled? */
   GLboolean StippleFlag;	/**< GL_LINE_STIPPLE enabled? */
   GLushort StipplePattern;	/**< Stipple pattern */
   GLint StippleFactor;		/**< Stipple repeat factor */
   GLfloat Width;		/**< Line width */
   GLfloat _Width;		/**< Clamped Line width */
};

struct Mgl_list_attrib {
   GLuint ListBase;
};

struct Mgl_multisample_attrib {
   GLboolean Enabled;
   GLboolean SampleAlphaToCoverage;
   GLboolean SampleAlphaToOne;
   GLboolean SampleCoverage;
   GLfloat SampleCoverageValue;
   GLboolean SampleCoverageInvert;
};

/**
 * Data structure for color tables
 */
struct Mgl_color_table {
   GLenum Format;         /**< GL_ALPHA, GL_RGB, GL_RGB, etc */
   GLenum IntFormat;
   GLuint Size;           /**< number of entries (rows) in table */
   GLvoid *Table;         /**< points to data of <Type> */
   GLenum Type;           /**< GL_UNSIGNED_BYTE or GL_FLOAT */
   GLubyte RedSize;
   GLubyte GreenSize;
   GLubyte BlueSize;
   GLubyte AlphaSize;
   GLubyte LuminanceSize;
   GLubyte IntensitySize;
};

/**
 * Pixel attributes.
 */
struct Mgl_pixel_attrib {
   GLenum ReadBuffer;		/**< source buffer for glReadPixels()/glCopyPixels() */
   GLubyte _ReadSrcMask;	/**< Not really a mask, but like _DrawDestMask
				  *
				  * May be: FRONT_LEFT_BIT, BACK_LEFT_BIT,
				  * FRONT_RIGHT_BIT or BACK_RIGHT_BIT. */
   GLfloat RedBias, RedScale;
   GLfloat GreenBias, GreenScale;
   GLfloat BlueBias, BlueScale;
   GLfloat AlphaBias, AlphaScale;
   GLfloat DepthBias, DepthScale;
   GLint IndexShift, IndexOffset;
   GLboolean MapColorFlag;
   GLboolean MapStencilFlag;
   GLfloat ZoomX, ZoomY;
   /* XXX move these out of gl_pixel_attrib */
   GLint MapStoSsize;		/**< Size of each pixel map */
   GLint MapItoIsize;
   GLint MapItoRsize;
   GLint MapItoGsize;
   GLint MapItoBsize;
   GLint MapItoAsize;
   GLint MapRtoRsize;
   GLint MapGtoGsize;
   GLint MapBtoBsize;
   GLint MapAtoAsize;
   GLint MapStoS[MMAX_PIXEL_MAP_TABLE];	/**< Pixel map tables */
   GLint MapItoI[MMAX_PIXEL_MAP_TABLE];
   GLfloat MapItoR[MMAX_PIXEL_MAP_TABLE];
   GLfloat MapItoG[MMAX_PIXEL_MAP_TABLE];
   GLfloat MapItoB[MMAX_PIXEL_MAP_TABLE];
   GLfloat MapItoA[MMAX_PIXEL_MAP_TABLE];
   GLubyte MapItoR8[MMAX_PIXEL_MAP_TABLE];  /**< converted to 8-bit color */
   GLubyte MapItoG8[MMAX_PIXEL_MAP_TABLE];
   GLubyte MapItoB8[MMAX_PIXEL_MAP_TABLE];
   GLubyte MapItoA8[MMAX_PIXEL_MAP_TABLE];
   GLfloat MapRtoR[MMAX_PIXEL_MAP_TABLE];
   GLfloat MapGtoG[MMAX_PIXEL_MAP_TABLE];
   GLfloat MapBtoB[MMAX_PIXEL_MAP_TABLE];
   GLfloat MapAtoA[MMAX_PIXEL_MAP_TABLE];
   /** GL_EXT_histogram */
   GLboolean HistogramEnabled;
   GLboolean MinMaxEnabled;
   /** GL_SGIS_pixel_texture */
   GLboolean PixelTextureEnabled;
   GLenum FragmentRgbSource;
   GLenum FragmentAlphaSource;
   /** GL_SGI_color_matrix */
   GLfloat PostColorMatrixScale[4];  /**< RGBA */
   GLfloat PostColorMatrixBias[4];   /**< RGBA */
   /** GL_SGI_color_table */
   GLfloat ColorTableScale[4];
   GLfloat ColorTableBias[4];
   GLboolean ColorTableEnabled;
   GLfloat PCCTscale[4];
   GLfloat PCCTbias[4];
   GLboolean PostConvolutionColorTableEnabled;
   GLfloat PCMCTscale[4];
   GLfloat PCMCTbias[4];
   GLboolean PostColorMatrixColorTableEnabled;
   /** GL_SGI_texture_color_table */
   GLfloat TextureColorTableScale[4];
   GLfloat TextureColorTableBias[4];
   /** Convolution */
   GLboolean Convolution1DEnabled;
   GLboolean Convolution2DEnabled;
   GLboolean Separable2DEnabled;
   GLfloat ConvolutionBorderColor[3][4];
   GLenum ConvolutionBorderMode[3];
   GLfloat ConvolutionFilterScale[3][4];
   GLfloat ConvolutionFilterBias[3][4];
   GLfloat PostConvolutionScale[4];  /**< RGBA */
   GLfloat PostConvolutionBias[4];   /**< RGBA */
};

/**
 * Point attributes.
 */
struct Mgl_point_attrib {
   GLboolean SmoothFlag;	/**< True if GL_POINT_SMOOTH is enabled */
   GLfloat Size;		/**< User-specified point size */
   GLfloat _Size;		/**< Size clamped to Const.Min/MaxPointSize */
   GLfloat Params[3];		/**< GL_EXT_point_parameters */
   GLfloat MinSize, MaxSize;	/**< GL_EXT_point_parameters */
   GLfloat Threshold;		/**< GL_EXT_point_parameters */
   GLboolean _Attenuated;	/**< True if Params != [1, 0, 0] */
   GLboolean PointSprite;	/**< GL_NV_point_sprite / GL_NV_point_sprite */
   GLboolean CoordReplace[MMAX_TEXTURE_UNITS]; /**< GL_NV_point_sprite / GL_NV_point_sprite */
   GLenum SpriteRMode;		/**< GL_NV_point_sprite (only!) */
};

/**
 * Polygon attributes.
 */
struct Mgl_polygon_attrib {
   GLenum FrontFace;		/**< Either GL_CW or GL_CCW */
   GLenum FrontMode;		/**< Either GL_POINT, GL_LINE or GL_FILL */
   GLenum BackMode;		/**< Either GL_POINT, GL_LINE or GL_FILL */
   GLboolean _FrontBit;		/**< 0=GL_CCW, 1=GL_CW */
   GLboolean CullFlag;		/**< Culling on/off flag */
   GLboolean SmoothFlag;	/**< True if GL_POLYGON_SMOOTH is enabled */
   GLboolean StippleFlag;	/**< True if GL_POLYGON_STIPPLE is enabled */
   GLenum CullFaceMode;		/**< Culling mode GL_FRONT or GL_BACK */
   GLfloat OffsetFactor;	/**< Polygon offset factor, from user */
   GLfloat OffsetUnits;		/**< Polygon offset units, from user */
   GLboolean OffsetPoint;	/**< Offset in GL_POINT mode */
   GLboolean OffsetLine;	/**< Offset in GL_LINE mode */
   GLboolean OffsetFill;	/**< Offset in GL_FILL mode */
};

/**
 * Scissor attributes.
 */
struct Mgl_scissor_attrib {
   GLboolean Enabled;		/**< Scissor test enabled? */
   GLint X, Y;			/**< Lower left corner of box */
   GLsizei Width, Height;	/**< Size of box */
};

/**
 * Stencil attributes.
 */
struct Mgl_stencil_attrib {
	GLboolean Enabled;		/**< Enabled flag */
	GLboolean TestTwoSide;	/**< GL_EXT_stencil_two_side */
	GLubyte ActiveFace;		/**< GL_EXT_stencil_two_side (0 or 1) */
	GLenum Function[2];		/**< Stencil function */
	GLenum FailFunc[2];		/**< Fail function */
	GLenum ZPassFunc[2];		/**< Depth buffer pass function */
	GLenum ZFailFunc[2];		/**< Depth buffer fail function */
	MGLstencil Ref[2];		/**< Reference value */
	MGLstencil ValueMask[2];	/**< Value mask */
	MGLstencil WriteMask[2];	/**< Write mask */
	MGLstencil Clear;		/**< Clear value */
};

/**
 * Texture object record
 */
struct Mgl_texture_object {
	_glthread_Mutex Mutex;	/**< for thread safety */
	GLint RefCount;		/**< reference count */
	GLuint Name;			/**< an unsigned integer */
	GLenum Target;               /**< GL_TEXTURE_1D, GL_TEXTURE_2D, etc. */
	GLfloat Priority;		/**< in [0,1] */
	GLfloat BorderColor[4];	/**< unclamped */
	MGLchan _BorderChan[4];	/**< clamped, as GLchan */
	/** \name Wrap modes
	* Are GL_CLAMP, REPEAT, GL_CLAMP_TO_EDGE, and GL_CLAMP_TO_BORDER_ARB. */
	/*@{*/
	GLenum WrapS;
	GLenum WrapT;
	GLenum WrapR;
	/*@}*/
	GLenum MinFilter;		/**< minification filter */
	GLenum MagFilter;		/**< magnification filter */
	GLfloat MinLod;		/**< min lambda, OpenGL 1.2 */
	GLfloat MaxLod;		/**< max lambda, OpenGL 1.2 */
	GLfloat LodBias;		/**< OpenGL 1.4 */
	GLint BaseLevel;		/**< min mipmap level, OpenGL 1.2 */
	GLint MaxLevel;		/**< max mipmap level, OpenGL 1.2 */
	GLfloat MaxAnisotropy;	/**< GL_EXT_texture_filter_anisotropic */
	GLboolean CompareFlag;	/**< GL_SGIX_shadow */
	GLenum CompareOperator;	/**< GL_SGIX_shadow */
	GLfloat ShadowAmbient;
	GLenum CompareMode;		/**< GL_ARB_shadow */
	GLenum CompareFunc;		/**< GL_ARB_shadow */
	GLenum DepthMode;		/**< GL_ARB_depth_texture */
	GLint _MaxLevel;		/**< actual max mipmap level (q in the spec) */
	GLfloat _MaxLambda;		/**< = _MaxLevel - BaseLevel (q - b in spec) */
	GLboolean GenerateMipmap;    /**< GL_SGIS_generate_mipmap */
	GLboolean _IsPowerOfTwo;	/**< Are all image dimensions powers of two? */

	struct gl_texture_image *Image[MMAX_FACES][MMAX_TEXTURE_LEVELS];

	/** GL_EXT_paletted_texture */
	struct Mgl_color_table Palette;

	GLboolean Complete;			/**< Is texture object complete? */
	struct gl_texture_object *Next;	/**< Next in linked list */

	/**
	* \name For device driver
	*/
	/*@{*/
	void *DriverData;	/**< Arbitrary device driver data */
	/*@}*/
};

/**
 * Texture combine environment state.
 * 
 * \todo
 * If GL_NV_texture_env_combine4 is ever supported, the arrays in this
 * structure will need to be expanded for 4 elements.
 */
struct Mgl_tex_env_combine_state {
   GLenum ModeRGB;       /**< GL_REPLACE, GL_DECAL, GL_ADD, etc. */
   GLenum ModeA;         /**< GL_REPLACE, GL_DECAL, GL_ADD, etc. */
   GLenum SourceRGB[3];  /**< GL_PRIMARY_COLOR, GL_TEXTURE, etc. */
   GLenum SourceA[3];    /**< GL_PRIMARY_COLOR, GL_TEXTURE, etc. */
   GLenum OperandRGB[3]; /**< SRC_COLOR, ONE_MINUS_SRC_COLOR, etc */
   GLenum OperandA[3];   /**< SRC_ALPHA, ONE_MINUS_SRC_ALPHA, etc */
   GLuint ScaleShiftRGB; /**< 0, 1 or 2 */
   GLuint ScaleShiftA;   /**< 0, 1 or 2 */
   GLuint _NumArgsRGB;   /**< Number of inputs used for the combine mode. */
   GLuint _NumArgsA;     /**< Number of inputs used for the combine mode. */
};

/**
 * Texture unit record 
 */
struct Mgl_texture_unit {
   GLuint Enabled;              /**< bitmask of TEXTURE_*_BIT flags */
   GLuint _ReallyEnabled;       /**< 0 or exactly one of TEXTURE_*_BIT flags */

   GLenum EnvMode;              /**< GL_MODULATE, GL_DECAL, GL_BLEND, etc. */
   GLfloat EnvColor[4];
   GLuint TexGenEnabled;	/**< Bitwise-OR of [STRQ]_BIT values */
   /** \name Tex coord generation mode
    * Either GL_OBJECT_LINEAR, GL_EYE_LINEAR or GL_SPHERE_MAP. */
   /*@{*/
   GLenum GenModeS;		
   GLenum GenModeT;
   GLenum GenModeR;
   GLenum GenModeQ;
   /*@}*/
   GLuint _GenBitS;
   GLuint _GenBitT;
   GLuint _GenBitR;
   GLuint _GenBitQ;
   GLuint _GenFlags;		/**< bitwise or of GenBit[STRQ] */
   GLfloat ObjectPlaneS[4];
   GLfloat ObjectPlaneT[4];
   GLfloat ObjectPlaneR[4];
   GLfloat ObjectPlaneQ[4];
   GLfloat EyePlaneS[4];
   GLfloat EyePlaneT[4];
   GLfloat EyePlaneR[4];
   GLfloat EyePlaneQ[4];
   GLfloat LodBias;		/**< for biasing mipmap levels */

   /** 
    * \name GL_EXT_texture_env_combine 
    */
   struct Mgl_tex_env_combine_state Combine;

   /**
    * Derived state based on \c EnvMode and the \c BaseFormat of the
    * currently enabled texture.
    */
   struct Mgl_tex_env_combine_state _EnvMode;

   /**
    * Currently enabled combiner state.  This will point to either
    * \c Combine or \c _EnvMode.
    */
   struct gl_tex_env_combine_state *_CurrentCombine;

   struct gl_texture_object *Current1D;
   struct gl_texture_object *Current2D;
   struct gl_texture_object *Current3D;
   struct gl_texture_object *CurrentCubeMap; /**< GL_ARB_texture_cube_map */
   struct gl_texture_object *CurrentRect;    /**< GL_NV_texture_rectangle */

   struct gl_texture_object *_Current; /**< Points to really enabled tex obj */

   struct Mgl_texture_object Saved1D;  /**< only used by glPush/PopAttrib */
   struct Mgl_texture_object Saved2D;
   struct Mgl_texture_object Saved3D;
   struct Mgl_texture_object SavedCubeMap;
   struct Mgl_texture_object SavedRect;

   /* GL_SGI_texture_color_table */
   struct Mgl_color_table ColorTable;
   struct Mgl_color_table ProxyColorTable;
   GLboolean ColorTableEnabled;
};

/**
 * Texture attributes
 */
struct Mgl_texture_attrib {
   /**
    * name multitexture 
    */
   /**@{*/
   GLuint CurrentUnit;	        /**< Active texture unit */
   GLuint _EnabledUnits;        /**< one bit set for each really-enabled unit */
   GLuint _EnabledCoordUnits;   /**< one bit per enabled coordinate unit */
   GLuint _GenFlags;            /**< for texgen */
   GLuint _TexGenEnabled;
   GLuint _TexMatEnabled;
   /**@}*/

   struct Mgl_texture_unit Unit[MMAX_TEXTURE_UNITS];

   struct gl_texture_object *Proxy1D;
   struct gl_texture_object *Proxy2D;
   struct gl_texture_object *Proxy3D;
   struct gl_texture_object *ProxyCubeMap;
   struct gl_texture_object *ProxyRect;

   /** GL_EXT_shared_texture_palette */
   GLboolean SharedPalette;
   struct Mgl_color_table Palette;
};

/**
 * Histogram attributes.
 */
struct Mgl_histogram_attrib {
   GLuint Width;				/**< number of table entries */
   GLint Format;				/**< GL_ALPHA, GL_RGB, etc */
   GLuint Count[MHISTOGRAM_TABLE_SIZE][4];	/**< the histogram */
   GLboolean Sink;				/**< terminate image transfer? */
   GLubyte RedSize;				/**< Bits per counter */
   GLubyte GreenSize;
   GLubyte BlueSize;
   GLubyte AlphaSize;
   GLubyte LuminanceSize;
};

struct Mgl_minmax_attrib {
   GLenum Format;
   GLboolean Sink;
   GLfloat Min[4], Max[4];   /**< RGBA */
};


struct Mgl_convolution_attrib {
   GLenum Format;
   GLenum InternalFormat;
   GLuint Width;
   GLuint Height;
   GLfloat Filter[MMAX_CONVOLUTION_WIDTH * MMAX_CONVOLUTION_HEIGHT * 4];
};

/**
 * Client vertex array attributes
 */
struct Mgl_client_array {
   GLint Size;                  /**< components per element (1,2,3,4) */
   GLenum Type;                 /**< datatype: GL_FLOAT, GL_INT, etc */
   GLsizei Stride;		/**< user-specified stride */
   GLsizei StrideB;		/**< actual stride in bytes */
   const GLubyte *Ptr;          /**< Points to array data */
   GLuint Enabled;		/**< one of the _NEW_ARRAY_ bits */
   GLboolean Normalized;        /**< GL_ARB_vertex_program */

   /**< GL_ARB_vertex_buffer_object */
   struct Mgl_buffer_object *BufferObj;
   GLuint _MaxElement;

   GLuint Flags;
};

/**
 * Array attributes.
 */
struct Mgl_array_attrib {
   struct Mgl_client_array Vertex;	     /**< client data descriptors */
   struct Mgl_client_array Normal;
   struct Mgl_client_array Color;
   struct Mgl_client_array SecondaryColor;
   struct Mgl_client_array FogCoord;
   struct Mgl_client_array Index;
   struct Mgl_client_array TexCoord[MMAX_TEXTURE_COORD_UNITS];
   struct Mgl_client_array EdgeFlag;

   struct Mgl_client_array VertexAttrib[MVERT_ATTRIB_MAX];  /**< GL_NV_vertex_program */

   GLint ActiveTexture;		/**< Client Active Texture */
   GLuint LockFirst;            /**< GL_EXT_compiled_vertex_array */
   GLuint LockCount;            /**< GL_EXT_compiled_vertex_array */

   GLuint _Enabled;		/**< _NEW_ARRAY_* - bit set if array enabled */
   GLuint NewState;		/**< _NEW_ARRAY_* */

#if FEATURE_ARB_vertex_buffer_object
   struct Mgl_buffer_object *NullBufferObj;
   struct Mgl_buffer_object *ArrayBufferObj;
   struct Mgl_buffer_object *ElementArrayBufferObj;
#endif
   GLuint _MaxElement;          /* Min of all enabled array's maxes */
};

/**
 * Client pixel packing/unpacking attributes
 */
struct Mgl_pixelstore_attrib {
   GLint Alignment;
   GLint RowLength;
   GLint SkipPixels;
   GLint SkipRows;
   GLint ImageHeight;     /**< for GL_EXT_texture3D */
   GLint SkipImages;      /**< for GL_EXT_texture3D */
   GLboolean SwapBytes;
   GLboolean LsbFirst;
   GLboolean ClientStorage; /**< GL_APPLE_client_storage */
   GLboolean Invert;        /**< GL_MESA_pack_invert */
   struct Mgl_buffer_object *BufferObj; /**< GL_ARB_pixel_buffer_object */
};

/**
 * 1-D Evaluator control points
 */
struct Mgl_1d_map
{
   GLuint Order;	/**< Number of control points */
   GLfloat u1, u2, du;	/**< u1, u2, 1.0/(u2-u1) */
   GLfloat *Points;	/**< Points to contiguous control points */
};

/**
 * 2-D Evaluator control points
 */
struct Mgl_2d_map
{
   GLuint Uorder;		/**< Number of control points in U dimension */
   GLuint Vorder;		/**< Number of control points in V dimension */
   GLfloat u1, u2, du;
   GLfloat v1, v2, dv;
   GLfloat *Points;		/**< Points to contiguous control points */
};


/**
 * All evaluator control points
 */
struct Mgl_evaluators
{
   /** 
    * \name 1-D maps
    */
   /*@{*/
   struct Mgl_1d_map Map1Vertex3;
   struct Mgl_1d_map Map1Vertex4;
   struct Mgl_1d_map Map1Index;
   struct Mgl_1d_map Map1Color4;
   struct Mgl_1d_map Map1Normal;
   struct Mgl_1d_map Map1Texture1;
   struct Mgl_1d_map Map1Texture2;
   struct Mgl_1d_map Map1Texture3;
   struct Mgl_1d_map Map1Texture4;
   struct Mgl_1d_map Map1Attrib[16];  /**< GL_NV_vertex_program */
   /*@}*/

   /** 
    * \name 2-D maps 
    */
   /*@{*/
   struct Mgl_2d_map Map2Vertex3;
   struct Mgl_2d_map Map2Vertex4;
   struct Mgl_2d_map Map2Index;
   struct Mgl_2d_map Map2Color4;
   struct Mgl_2d_map Map2Normal;
   struct Mgl_2d_map Map2Texture1;
   struct Mgl_2d_map Map2Texture2;
   struct Mgl_2d_map Map2Texture3;
   struct Mgl_2d_map Map2Texture4;
   struct Mgl_2d_map Map2Attrib[16];  /**< GL_NV_vertex_program */
   /*@}*/
};

struct Mgl_feedback {
   GLenum Type;
   GLuint _Mask;		/* FB_* bits */
   GLfloat *Buffer;
   GLuint BufferSize;
   GLuint Count;
};

/**
 * Selection attributes.
 */
struct Mgl_selection {
   GLuint *Buffer;	/**< selection buffer */
   GLuint BufferSize;	/**< size of the selection buffer */
   GLuint BufferCount;	/**< number of values in the selection buffer */
   GLuint Hits;		/**< number of records in the selection buffer */
   GLuint NameStackDepth; /**< name stack depth */
   GLuint NameStack[MMAX_NAME_STACK_DEPTH]; /**< name stack */
   GLboolean HitFlag;	/**< hit flag */
   GLfloat HitMinZ;	/**< minimum hit depth */
   GLfloat HitMaxZ;	/**< maximum hit depth */
};

/**
 * State common to vertex and fragment programs.
 */
struct Mprogram_state {
   GLint ErrorPos;                       /* GL_PROGRAM_ERROR_POSITION_NV */
   const char *ErrorString;              /* GL_PROGRAM_ERROR_STRING_NV */
};

/**
 * State vars for GL_NV_vertex_program
 */
struct Mvertex_program_state
{
   GLboolean Enabled;                  /**< GL_VERTEX_PROGRAM_NV */
   GLboolean _Enabled;                 /**< Really enabled? */
   GLboolean PointSizeEnabled;         /**< GL_VERTEX_PROGRAM_POINT_SIZE_NV */
   GLboolean TwoSideEnabled;           /**< GL_VERTEX_PROGRAM_TWO_SIDE_NV */
   struct vertex_program *Current;     /**< ptr to currently bound program */

   GLenum TrackMatrix[MMAX_NV_VERTEX_PROGRAM_PARAMS / 4];
   GLenum TrackMatrixTransform[MMAX_NV_VERTEX_PROGRAM_PARAMS / 4];

   GLfloat Parameters[MMAX_NV_VERTEX_PROGRAM_PARAMS][4]; /* Env params */
   /* Only used during program execution (may be moved someday): */
   GLfloat Temporaries[MMAX_NV_VERTEX_PROGRAM_TEMPS][4];
   GLfloat Inputs[MMAX_NV_VERTEX_PROGRAM_INPUTS][4];
   GLfloat Outputs[MMAX_NV_VERTEX_PROGRAM_OUTPUTS][4];
   GLint AddressReg[4];

#if FEATURE_MESA_program_debug
   GLprogramcallbackMESA Callback;
   GLvoid *CallbackData;
   GLboolean CallbackEnabled;
   GLuint CurrentPosition;
#endif
};

/**
 * Base class for any kind of program object
 */
struct Mprogram
{
   GLuint Id;
   GLubyte *String;    /* Null-terminated program text */
   GLenum Target;
   GLenum Format;      /* String encoding format */
   GLint RefCount;
   GLboolean Resident;
   GLfloat LocalParams[MMAX_PROGRAM_LOCAL_PARAMS][4];
   GLuint NumInstructions;  /* GL_ARB_vertex/fragment_program */
   GLuint NumTemporaries;
   GLuint NumParameters;
   GLuint NumAttributes;
   GLuint NumAddressRegs;
};

/** Fragment program object */
struct Mfragment_program
{
   struct Mprogram Base;   /**< base class */
   struct fp_instruction *Instructions;  /**< Compiled instructions */
   GLuint InputsRead;     /**< Bitmask of which input regs are read */
   GLuint OutputsWritten; /**< Bitmask of which output regs are written to */
   GLuint TexturesUsed[MMAX_TEXTURE_IMAGE_UNITS];  /**< TEXTURE_x_INDEX bitmask */
   GLuint NumAluInstructions; /**< GL_ARB_fragment_program */
   GLuint NumTexInstructions;
   GLuint NumTexIndirections;
   GLenum FogOption;
   struct program_parameter_list *Parameters; /**< array [NumParameters] */

#ifdef USE_TCC
   char c_str[4096];		/* experimental... */
   int c_strlen;
#endif
};

/**
 * NV_fragment_program runtime state
 */
struct Mfp_machine
{
   GLfloat Temporaries[MMAX_NV_FRAGMENT_PROGRAM_TEMPS][4];
   GLfloat Inputs[MMAX_NV_FRAGMENT_PROGRAM_INPUTS][4];
   GLfloat Outputs[MMAX_NV_FRAGMENT_PROGRAM_OUTPUTS][4];
   GLuint CondCodes[4];
};

/*
 * State for GL_ARB/NV_fragment_program
 */
struct Mfragment_program_state
{
   GLboolean Enabled;                    /* GL_VERTEX_PROGRAM_NV */
   GLboolean _Enabled;                   /* Really enabled? */
   struct Mfragment_program *Current;     /* ptr to currently bound program */
   struct Mfp_machine Machine;            /* machine state */
   GLfloat Parameters[MMAX_NV_FRAGMENT_PROGRAM_PARAMS][4]; /* Env params */

#if FEATURE_MESA_program_debug
   GLprogramcallbackMESA Callback;
   GLvoid *CallbackData;
   GLboolean CallbackEnabled;
   GLuint CurrentPosition;
#endif
};

/**
 * An entry in the hash table.  
 *
 * This struct is private to this file.
 */
struct MHashEntry {
   GLuint Key;             /**< the entry's key */
   void *Data;             /**< the entry's data */
   struct HashEntry *Next; /**< pointer to next entry */
};

/**
 * The hash table data structure.  
 *
 * This is an opaque types (it's not defined in hash.h file).
 */
struct M_mesa_HashTable {
   struct MHashEntry *Table[MTABLE_SIZE];  /**< the lookup table */
   GLuint MaxKey;                        /**< highest key inserted so far */
   _glthread_Mutex Mutex;                /**< mutual exclusion lock */
};

/*
 * State for GL_ARB_occlusion_query
 */
struct Mocclusion_state
{
   GLboolean Active;
   GLuint CurrentQueryObject;
   GLuint PassedCounter;
   struct _mesa_HashTable *QueryObjects;
};

struct Mgl_list_instruction {
   GLuint Size;
   void (*Execute)( MGLcontext *ctx, void *data );
   void (*Destroy)( MGLcontext *ctx, void *data );
   void (*Print)( MGLcontext *ctx, void *data );
};

struct Mgl_list_extensions {
   struct Mgl_list_instruction Opcode[MMAX_DLIST_EXT_OPCODES];
   GLuint NumOpcodes;
};

/**
 * Transform/Clip/Lighting interface
 *
 * Drivers present a reduced set of the functions possible in
 * glBegin()/glEnd() objects.  Core mesa provides translation stubs for the
 * remaining functions to map down to these entry points.
 *
 * These are the initial values to be installed into dispatch by
 * mesa.  If the T&L driver wants to modify the dispatch table
 * while installed, it must do so itself.  It would be possible for
 * the vertexformat to install it's own initial values for these
 * functions, but this way there is an obvious list of what is
 * expected of the driver.
 *
 * If the driver wants to hook in entry points other than those
 * listed, it must restore them to their original values in
 * the disable() callback, below.
 */
typedef struct {
   /**
    * \name Vertex
    */
   /*@{*/
   void (GLAPIENTRYP ArrayElement)( GLint ); /* NOTE */
   void (GLAPIENTRYP Color3f)( GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP Color3fv)( const GLfloat * );
   void (GLAPIENTRYP Color4f)( GLfloat, GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP Color4fv)( const GLfloat * );
   void (GLAPIENTRYP EdgeFlag)( GLboolean );
   void (GLAPIENTRYP EdgeFlagv)( const GLboolean * );
   void (GLAPIENTRYP EvalCoord1f)( GLfloat );          /* NOTE */
   void (GLAPIENTRYP EvalCoord1fv)( const GLfloat * ); /* NOTE */
   void (GLAPIENTRYP EvalCoord2f)( GLfloat, GLfloat ); /* NOTE */
   void (GLAPIENTRYP EvalCoord2fv)( const GLfloat * ); /* NOTE */
   void (GLAPIENTRYP EvalPoint1)( GLint );             /* NOTE */
   void (GLAPIENTRYP EvalPoint2)( GLint, GLint );      /* NOTE */
   void (GLAPIENTRYP FogCoordfEXT)( GLfloat );
   void (GLAPIENTRYP FogCoordfvEXT)( const GLfloat * );
   void (GLAPIENTRYP Indexf)( GLfloat );
   void (GLAPIENTRYP Indexfv)( const GLfloat * );
   void (GLAPIENTRYP Materialfv)( GLenum face, GLenum pname, const GLfloat * ); /* NOTE */
   void (GLAPIENTRYP MultiTexCoord1fARB)( GLenum, GLfloat );
   void (GLAPIENTRYP MultiTexCoord1fvARB)( GLenum, const GLfloat * );
   void (GLAPIENTRYP MultiTexCoord2fARB)( GLenum, GLfloat, GLfloat );
   void (GLAPIENTRYP MultiTexCoord2fvARB)( GLenum, const GLfloat * );
   void (GLAPIENTRYP MultiTexCoord3fARB)( GLenum, GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP MultiTexCoord3fvARB)( GLenum, const GLfloat * );
   void (GLAPIENTRYP MultiTexCoord4fARB)( GLenum, GLfloat, GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP MultiTexCoord4fvARB)( GLenum, const GLfloat * );
   void (GLAPIENTRYP Normal3f)( GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP Normal3fv)( const GLfloat * );
   void (GLAPIENTRYP SecondaryColor3fEXT)( GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP SecondaryColor3fvEXT)( const GLfloat * );
   void (GLAPIENTRYP TexCoord1f)( GLfloat );
   void (GLAPIENTRYP TexCoord1fv)( const GLfloat * );
   void (GLAPIENTRYP TexCoord2f)( GLfloat, GLfloat );
   void (GLAPIENTRYP TexCoord2fv)( const GLfloat * );
   void (GLAPIENTRYP TexCoord3f)( GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP TexCoord3fv)( const GLfloat * );
   void (GLAPIENTRYP TexCoord4f)( GLfloat, GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP TexCoord4fv)( const GLfloat * );
   void (GLAPIENTRYP Vertex2f)( GLfloat, GLfloat );
   void (GLAPIENTRYP Vertex2fv)( const GLfloat * );
   void (GLAPIENTRYP Vertex3f)( GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP Vertex3fv)( const GLfloat * );
   void (GLAPIENTRYP Vertex4f)( GLfloat, GLfloat, GLfloat, GLfloat );
   void (GLAPIENTRYP Vertex4fv)( const GLfloat * );
   void (GLAPIENTRYP CallList)( GLuint );	/* NOTE */
   void (GLAPIENTRYP CallLists)( GLsizei, GLenum, const GLvoid * );	/* NOTE */
   void (GLAPIENTRYP Begin)( GLenum );
   void (GLAPIENTRYP End)( void );
   void (GLAPIENTRYP VertexAttrib1fNV)( GLuint index, GLfloat x );
   void (GLAPIENTRYP VertexAttrib1fvNV)( GLuint index, const GLfloat *v );
   void (GLAPIENTRYP VertexAttrib2fNV)( GLuint index, GLfloat x, GLfloat y );
   void (GLAPIENTRYP VertexAttrib2fvNV)( GLuint index, const GLfloat *v );
   void (GLAPIENTRYP VertexAttrib3fNV)( GLuint index, GLfloat x, GLfloat y, GLfloat z );
   void (GLAPIENTRYP VertexAttrib3fvNV)( GLuint index, const GLfloat *v );
   void (GLAPIENTRYP VertexAttrib4fNV)( GLuint index, GLfloat x, GLfloat y, GLfloat z, GLfloat w );
   void (GLAPIENTRYP VertexAttrib4fvNV)( GLuint index, const GLfloat *v );
   /*@}*/

   /*
    */
   void (GLAPIENTRYP Rectf)( GLfloat, GLfloat, GLfloat, GLfloat );

   /**
    * \name Array
    */
   /*@{*/
   void (GLAPIENTRYP DrawArrays)( GLenum mode, GLint start, GLsizei count );
   void (GLAPIENTRYP DrawElements)( GLenum mode, GLsizei count, GLenum type,
			 const GLvoid *indices );
   void (GLAPIENTRYP DrawRangeElements)( GLenum mode, GLuint start,
			      GLuint end, GLsizei count,
			      GLenum type, const GLvoid *indices );
   /*@}*/

   /**
    * \name Eval
    *
    * If you don't support eval, fallback to the default vertex format
    * on receiving an eval call and use the pipeline mechanism to
    * provide partial T&L acceleration.
    *
    * Mesa will provide a set of helper functions to do eval within
    * accelerated vertex formats, eventually...
    */
   /*@{*/
   void (GLAPIENTRYP EvalMesh1)( GLenum mode, GLint i1, GLint i2 );
   void (GLAPIENTRYP EvalMesh2)( GLenum mode, GLint i1, GLint i2, GLint j1, GLint j2 );
   /*@}*/

} MGLvertexformat;

#define MNUM_VERTEX_FORMAT_ENTRIES (sizeof(MGLvertexformat) / sizeof(void *))

struct Mgl_tnl_module {
   /**
    * Vertex format to be lazily swapped into current dispatch.
    */
   const MGLvertexformat *Current;

   /**
    * \name Record of functions swapped out.  
    * On restore, only need to swap these functions back in.
    */
   /*@{*/
   void *Swapped[MNUM_VERTEX_FORMAT_ENTRIES][2];
   GLuint SwapCount;
   /*@}*/
};

/**
 * List of extensions.
 */
struct Mgl_extensions
{
   /**
    * \name Flags to quickly test if certain extensions are available.
    * 
    * Not every extension needs to have such a flag, but it's encouraged.
    */
   /*@{*/
   GLboolean dummy;  /* don't remove this! */
   GLboolean ARB_depth_texture;
   GLboolean ARB_fragment_program;
   GLboolean ARB_half_float_pixel;
   GLboolean ARB_imaging;
   GLboolean ARB_multisample;
   GLboolean ARB_multitexture;
   GLboolean ARB_occlusion_query;
   GLboolean ARB_point_sprite;
   GLboolean ARB_shadow;
   GLboolean ARB_texture_border_clamp;
   GLboolean ARB_texture_compression;
   GLboolean ARB_texture_cube_map;
   GLboolean ARB_texture_env_combine;
   GLboolean ARB_texture_env_crossbar;
   GLboolean ARB_texture_env_dot3;
   GLboolean ARB_texture_float;
   GLboolean ARB_texture_mirrored_repeat;
   GLboolean ARB_texture_non_power_of_two;
   GLboolean ARB_transpose_matrix;
   GLboolean ARB_vertex_buffer_object;
   GLboolean ARB_vertex_program;
   GLboolean ARB_window_pos;
   GLboolean EXT_abgr;
   GLboolean EXT_bgra;
   GLboolean EXT_blend_color;
   GLboolean EXT_blend_equation_separate;
   GLboolean EXT_blend_func_separate;
   GLboolean EXT_blend_logic_op;
   GLboolean EXT_blend_minmax;
   GLboolean EXT_blend_subtract;
   GLboolean EXT_clip_volume_hint;
   GLboolean EXT_cull_vertex;
   GLboolean EXT_convolution;
   GLboolean EXT_compiled_vertex_array;
   GLboolean EXT_copy_texture;
   GLboolean EXT_depth_bounds_test;
   GLboolean EXT_draw_range_elements;
   GLboolean EXT_fog_coord;
   GLboolean EXT_histogram;
   GLboolean EXT_multi_draw_arrays;
   GLboolean EXT_paletted_texture;
   GLboolean EXT_packed_pixels;
   GLboolean EXT_pixel_buffer_object;
   GLboolean EXT_point_parameters;
   GLboolean EXT_polygon_offset;
   GLboolean EXT_rescale_normal;
   GLboolean EXT_shadow_funcs;
   GLboolean EXT_secondary_color;
   GLboolean EXT_separate_specular_color;
   GLboolean EXT_shared_texture_palette;
   GLboolean EXT_stencil_wrap;
   GLboolean EXT_stencil_two_side;
   GLboolean EXT_subtexture;
   GLboolean EXT_texture;
   GLboolean EXT_texture_object;
   GLboolean EXT_texture3D;
   GLboolean EXT_texture_compression_s3tc;
   GLboolean EXT_texture_env_add;
   GLboolean EXT_texture_env_combine;
   GLboolean EXT_texture_env_dot3;
   GLboolean EXT_texture_filter_anisotropic;
   GLboolean EXT_texture_lod_bias;
   GLboolean EXT_texture_mirror_clamp;
   GLboolean EXT_vertex_array;
   GLboolean EXT_vertex_array_set;
   /* vendor extensions */
   GLboolean APPLE_client_storage;
   GLboolean APPLE_packed_pixels;
   GLboolean ATI_texture_mirror_once;
   GLboolean ATI_texture_env_combine3;
   GLboolean HP_occlusion_test;
   GLboolean IBM_rasterpos_clip;
   GLboolean IBM_multimode_draw_arrays;
   GLboolean MESA_pack_invert;
   GLboolean MESA_packed_depth_stencil;
   GLboolean MESA_program_debug;
   GLboolean MESA_resize_buffers;
   GLboolean MESA_ycbcr_texture;
   GLboolean NV_blend_square;
   GLboolean NV_fragment_program;
   GLboolean NV_light_max_exponent;
   GLboolean NV_point_sprite;
   GLboolean NV_texgen_reflection;
   GLboolean NV_texture_rectangle;
   GLboolean NV_vertex_program;
   GLboolean NV_vertex_program1_1;
   GLboolean SGI_color_matrix;
   GLboolean SGI_color_table;
   GLboolean SGI_texture_color_table;
   GLboolean SGIS_generate_mipmap;
   GLboolean SGIS_pixel_texture;
   GLboolean SGIS_texture_edge_clamp;
   GLboolean SGIS_texture_lod;
   GLboolean SGIX_depth_texture;
   GLboolean SGIX_pixel_texture;
   GLboolean SGIX_shadow;
   GLboolean SGIX_shadow_ambient; /* or GL_ARB_shadow_ambient */
   GLboolean TDFX_texture_compression_FXT1;
   GLboolean S3_s3tc;
   /*@}*/
   /* The extension string */
   const GLubyte *String;
};

/**
 * Matrix.
 */
typedef struct {
   GLfloat *m;		/**< matrix, 16-byte aligned */
   GLfloat *inv;	/**< optional inverse, 16-byte aligned */
   GLuint flags;        /**< possible values determined by (of \link
			   MatFlags MAT_FLAG_* flags\endlink) */
   enum GLmatrixtype type;
} MGLmatrix;

/**
 * A stack of matrices (projection, modelview, color, texture, etc).
 */
struct Mmatrix_stack
{
   MGLmatrix *Top;      /**< points into Stack */
   MGLmatrix *Stack;    /**< array [MaxDepth] of GLmatrix */
   GLuint Depth;       /**< 0 <= Depth < MaxDepth */
   GLuint MaxDepth;    /**< size of Stack[] array */
   GLuint DirtyFlag;   /**< _NEW_MODELVIEW or _NEW_PROJECTION, for example */
};

/**
 * Material shininess lookup table.
 */
struct Mgl_shine_tab {
   struct Mgl_shine_tab *next, *prev;
   GLfloat tab[MSHINE_TABLE_SIZE+1];
   GLfloat shininess;
   GLuint refcount;
};

/**
 * Viewport attributes.
 */
struct Mgl_viewport_attrib {
	GLint X, Y;			/**< position */
	GLsizei Width, Height;	/**< size */
	GLfloat Near, Far;		/**< Depth buffer range */
	MGLmatrix _WindowMap;		/**< Mapping transformation as a matrix. */
};

struct Mmesa_list_state {
   GLuint CallDepth;		/**< Current recursion calling depth */
   MNode *CurrentListPtr;	/**< Head of list being compiled */
   GLuint CurrentListNum;	/**< Number of the list being compiled */
   MNode *CurrentBlock;		/**< Pointer to current block of nodes */
   GLuint CurrentPos;		/**< Index into current block of nodes */
   MGLvertexformat ListVtxfmt;

   GLubyte ActiveAttribSize[MVERT_ATTRIB_MAX];
   GLfloat CurrentAttrib[MVERT_ATTRIB_MAX][4];
   
   GLubyte ActiveMaterialSize[MMAT_ATTRIB_MAX];
   GLfloat CurrentMaterial[MMAT_ATTRIB_MAX][4];

   GLubyte ActiveIndex;
   GLfloat CurrentIndex;
   
   GLubyte ActiveEdgeFlag;
   GLboolean CurrentEdgeFlag;
};

/*
** This file defines the interface between the GL core and the surrounding
** "operating system" that supports it (currently the GLX or WGL extensions).
**
** Members (data and function pointers) are documented as imported or
** exported according to how they are used by the core rendering functions.
** Imported members are initialized by the "operating system" and used by
** the core functions.  Exported members are initialized by the core functions
** and used by the "operating system".
*/

/*
** Mode and limit information for a context.  This information is
** kept around in the context so that values can be used during
** command execution, and for returning information about the
** context to the application.
*/
typedef struct M__GLcontextModesRec {
    struct M__GLcontextModesRec * next;

    GLboolean rgbMode;
    GLboolean floatMode;
    GLboolean colorIndexMode;
    GLuint doubleBufferMode;
    GLuint stereoMode;

    GLboolean haveAccumBuffer;
    GLboolean haveDepthBuffer;
    GLboolean haveStencilBuffer;

    GLint redBits, greenBits, blueBits, alphaBits;	/* bits per comp */
    GLuint redMask, greenMask, blueMask, alphaMask;
    GLint rgbBits;		/* total bits for rgb */
    GLint indexBits;		/* total bits for colorindex */

    GLint accumRedBits, accumGreenBits, accumBlueBits, accumAlphaBits;
    GLint depthBits;
    GLint stencilBits;

    GLint numAuxBuffers;

    GLint level;

    GLint pixmapMode;

    /* GLX */
    GLint visualID;
    GLint visualType;     /**< One of the GLX X visual types. (i.e., 
			   * \c GLX_TRUE_COLOR, etc.)
			   */

    /* EXT_visual_rating / GLX 1.2 */
    GLint visualRating;

    /* EXT_visual_info / GLX 1.2 */
    GLint transparentPixel;
				/*    colors are floats scaled to ints */
    GLint transparentRed, transparentGreen, transparentBlue, transparentAlpha;
    GLint transparentIndex;

    /* ARB_multisample / SGIS_multisample */
    GLint sampleBuffers;
    GLint samples;

    /* SGIX_fbconfig / GLX 1.3 */
    GLint drawableType;
    GLint renderType;
    GLint xRenderable;
    GLint fbconfigID;

    /* SGIX_pbuffer / GLX 1.3 */
    GLint maxPbufferWidth;
    GLint maxPbufferHeight;
    GLint maxPbufferPixels;
    GLint optimalPbufferWidth;   /* Only for SGIX_pbuffer. */
    GLint optimalPbufferHeight;  /* Only for SGIX_pbuffer. */

    /* SGIX_visual_select_group */
    GLint visualSelectGroup;

    /* OML_swap_method */
    GLint swapMethod;

    GLint screen;
} __GLcontextModes;

/**
 * Device driver function table.
 * Core Mesa uses these function pointers to call into device drivers.
 * Most of these functions directly correspond to OpenGL state commands.
 * Core Mesa will call these functions after error checking has been done
 * so that the drivers don't have to worry about error testing.
 *
 * Vertex transformation/clipping/lighting is patched into the T&L module.
 * Rasterization functions are patched into the swrast module.
 *
 * Note: when new functions are added here, the drivers/common/driverfuncs.c
 * file should be updated too!!!
 */
struct Mdd_function_table {
   /**
    * Return a string as needed by glGetString().
    *
    * Only the GL_RENDERER token must be implemented.  Otherwise, NULL can be
    * returned.
    */
   const GLubyte * (*GetString)( MGLcontext *ctx, GLenum name );

   /**
    * Notify the driver after Mesa has made some internal state changes.  
    *
    * This is in addition to any state change callbacks Mesa may already have
    * made.
    */
   void (*UpdateState)( MGLcontext *ctx, GLuint new_state );

   /**
    * Get the width and height of the named buffer/window.
    *
    * Mesa uses this to determine when the driver's window size has changed.
    */
   void (*GetBufferSize)( MGLframebuffer *buffer,
                          GLuint *width, GLuint *height );

   /**
    * Resize the driver's depth/stencil/accum/back buffers to match the
    * size given in the GLframebuffer struct.  
    *
    * This is typically called when Mesa detects that a window size has changed.
    */
   void (*ResizeBuffers)( MGLframebuffer *buffer );

   /**
    * Called whenever an error is generated.  
    *
    * __GLcontextRec::ErrorValue contains the error value.
    */
   void (*Error)( MGLcontext *ctx );

   /**
    * This is called whenever glFinish() is called.
    */
   void (*Finish)( MGLcontext *ctx );

   /**
    * This is called whenever glFlush() is called.
    */
   void (*Flush)( MGLcontext *ctx );

   /**
    * Clear the color/depth/stencil/accum buffer(s).
    *
    * \param mask a bitmask of the DD_*_BIT values defined above that indicates
    * which buffers need to be cleared.
    * \param all if true then clear the whole buffer, else clear only the
    * region defined by <tt>(x, y, width, height)</tt>.
    * 
    * This function must obey the glColorMask(), glIndexMask() and
    * glStencilMask() settings!
    * Software Mesa can do masked clears if the device driver can't.
    */
   void (*Clear)( MGLcontext *ctx, GLbitfield mask, GLboolean all,
		  GLint x, GLint y, GLint width, GLint height );


   /**
    * \name For hardware accumulation buffer
    */
   /*@{*/
   /**
    * Execute glAccum command within the given scissor region.
    */
   void (*Accum)( MGLcontext *ctx, GLenum op, GLfloat value,
		  GLint xpos, GLint ypos, GLint width, GLint height );
   /*@}*/


   /**
    * \name glDraw(), glRead(), glCopyPixels() and glBitmap() functions
    */
   /*@{*/

   /**
    * This is called by glDrawPixels().
    *
    * \p unpack describes how to unpack the source image data.
    */
   void (*DrawPixels)( MGLcontext *ctx,
		       GLint x, GLint y, GLsizei width, GLsizei height,
		       GLenum format, GLenum type,
		       const struct gl_pixelstore_attrib *unpack,
		       const GLvoid *pixels );

   /**
    * Called by glReadPixels().
    */
   void (*ReadPixels)( MGLcontext *ctx,
		       GLint x, GLint y, GLsizei width, GLsizei height,
		       GLenum format, GLenum type,
		       const struct gl_pixelstore_attrib *unpack,
		       GLvoid *dest );

   /**
    * Do a glCopyPixels().  
    *
    * This function must respect all rasterization state, glPixelTransfer(),
    * glPixelZoom(), etc.
    */
   void (*CopyPixels)( MGLcontext *ctx,
                            GLint srcx, GLint srcy,
                            GLsizei width, GLsizei height,
                            GLint dstx, GLint dsty, GLenum type );

   /**
    * This is called by glBitmap().  
    *
    * Works the same as dd_function_table::DrawPixels, above.
    */
   void (*Bitmap)( MGLcontext *ctx,
		   GLint x, GLint y, GLsizei width, GLsizei height,
		   const struct gl_pixelstore_attrib *unpack,
		   const GLubyte *bitmap );
   /*@}*/

   
   /**
    * \name Texture image functions
    */
   /*@{*/

   /**
    * Choose texture format.
    * 
    * This is called by the \c _mesa_store_tex[sub]image[123]d() fallback
    * functions.  The driver should examine \p internalFormat and return a
    * pointer to an appropriate gl_texture_format.
    */
   const struct gl_texture_format *(*ChooseTextureFormat)( MGLcontext *ctx,
                      GLint internalFormat, GLenum srcFormat, GLenum srcType );

   /**
    * Called by glTexImage1D().
    * 
    * \param target user specified.
    * \param format user specified.
    * \param type user specified.
    * \param pixels user specified.
    * \param packing indicates the image packing of pixels.
    * \param texObj is the target texture object.
    * \param texImage is the target texture image.  It will have the texture \p
    * width, \p height, \p depth, \p border and \p internalFormat information.
    * 
    * \p retainInternalCopy is returned by this function and indicates whether
    * core Mesa should keep an internal copy of the texture image.
    *
    * Drivers should call a fallback routine from texstore.c if needed.
    */
   void (*TexImage1D)( MGLcontext *ctx, GLenum target, GLint level,
                       GLint internalFormat,
                       GLint width, GLint border,
                       GLenum format, GLenum type, const GLvoid *pixels,
                       const struct gl_pixelstore_attrib *packing,
                       struct gl_texture_object *texObj,
                       struct gl_texture_image *texImage );

   /**
    * Called by glTexImage2D().
    * 
    * \sa dd_function_table::TexImage1D.
    */
   void (*TexImage2D)( MGLcontext *ctx, GLenum target, GLint level,
                       GLint internalFormat,
                       GLint width, GLint height, GLint border,
                       GLenum format, GLenum type, const GLvoid *pixels,
                       const struct gl_pixelstore_attrib *packing,
                       struct gl_texture_object *texObj,
                       struct gl_texture_image *texImage );
   
   /**
    * Called by glTexImage3D().
    * 
    * \sa dd_function_table::TexImage1D.
    */
   void (*TexImage3D)( MGLcontext *ctx, GLenum target, GLint level,
                       GLint internalFormat,
                       GLint width, GLint height, GLint depth, GLint border,
                       GLenum format, GLenum type, const GLvoid *pixels,
                       const struct gl_pixelstore_attrib *packing,
                       struct gl_texture_object *texObj,
                       struct gl_texture_image *texImage );

   /**
    * Called by glTexSubImage1D().
    *
    * \param target user specified.
    * \param level user specified.
    * \param xoffset user specified.
    * \param yoffset user specified.
    * \param zoffset user specified.
    * \param width user specified.
    * \param height user specified.
    * \param depth user specified.
    * \param format user specified.
    * \param type user specified.
    * \param pixels user specified.
    * \param packing indicates the image packing of pixels.
    * \param texObj is the target texture object.
    * \param texImage is the target texture image.  It will have the texture \p
    * width, \p height, \p border and \p internalFormat information.
    *
    * The driver should use a fallback routine from texstore.c if needed.
    */
   void (*TexSubImage1D)( MGLcontext *ctx, GLenum target, GLint level,
                          GLint xoffset, GLsizei width,
                          GLenum format, GLenum type,
                          const GLvoid *pixels,
                          const struct gl_pixelstore_attrib *packing,
                          struct gl_texture_object *texObj,
                          struct gl_texture_image *texImage );
   
   /**
    * Called by glTexSubImage2D().
    *
    * \sa dd_function_table::TexSubImage1D.
    */
   void (*TexSubImage2D)( MGLcontext *ctx, GLenum target, GLint level,
                          GLint xoffset, GLint yoffset,
                          GLsizei width, GLsizei height,
                          GLenum format, GLenum type,
                          const GLvoid *pixels,
                          const struct gl_pixelstore_attrib *packing,
                          struct gl_texture_object *texObj,
                          struct gl_texture_image *texImage );
   
   /**
    * Called by glTexSubImage3D().
    *
    * \sa dd_function_table::TexSubImage1D.
    */
   void (*TexSubImage3D)( MGLcontext *ctx, GLenum target, GLint level,
                          GLint xoffset, GLint yoffset, GLint zoffset,
                          GLsizei width, GLsizei height, GLint depth,
                          GLenum format, GLenum type,
                          const GLvoid *pixels,
                          const struct gl_pixelstore_attrib *packing,
                          struct gl_texture_object *texObj,
                          struct gl_texture_image *texImage );

   /**
    * Called by glCopyTexImage1D().
    * 
    * Drivers should use a fallback routine from texstore.c if needed.
    */
   void (*CopyTexImage1D)( MGLcontext *ctx, GLenum target, GLint level,
                           GLenum internalFormat, GLint x, GLint y,
                           GLsizei width, GLint border );

   /**
    * Called by glCopyTexImage2D().
    * 
    * Drivers should use a fallback routine from texstore.c if needed.
    */
   void (*CopyTexImage2D)( MGLcontext *ctx, GLenum target, GLint level,
                           GLenum internalFormat, GLint x, GLint y,
                           GLsizei width, GLsizei height, GLint border );

   /**
    * Called by glCopyTexSubImage1D().
    * 
    * Drivers should use a fallback routine from texstore.c if needed.
    */
   void (*CopyTexSubImage1D)( MGLcontext *ctx, GLenum target, GLint level,
                              GLint xoffset,
                              GLint x, GLint y, GLsizei width );
   /**
    * Called by glCopyTexSubImage2D().
    * 
    * Drivers should use a fallback routine from texstore.c if needed.
    */
   void (*CopyTexSubImage2D)( MGLcontext *ctx, GLenum target, GLint level,
                              GLint xoffset, GLint yoffset,
                              GLint x, GLint y,
                              GLsizei width, GLsizei height );
   /**
    * Called by glCopyTexSubImage3D().
    * 
    * Drivers should use a fallback routine from texstore.c if needed.
    */
   void (*CopyTexSubImage3D)( MGLcontext *ctx, GLenum target, GLint level,
                              GLint xoffset, GLint yoffset, GLint zoffset,
                              GLint x, GLint y,
                              GLsizei width, GLsizei height );

   /**
    * Called by glTexImage[123]D when user specifies a proxy texture
    * target.  
    *
    * \return GL_TRUE if the proxy test passes, or GL_FALSE if the test fails.
    */
   GLboolean (*TestProxyTexImage)(MGLcontext *ctx, GLenum target,
                                  GLint level, GLint internalFormat,
                                  GLenum format, GLenum type,
                                  GLint width, GLint height,
                                  GLint depth, GLint border);
   /*@}*/

   
   /**
    * \name Compressed texture functions
    */
   /*@{*/

   /**
    * Called by glCompressedTexImage1D().
    *
    * \param target user specified.
    * \param format user specified.
    * \param type user specified.
    * \param pixels user specified.
    * \param packing indicates the image packing of pixels.
    * \param texObj is the target texture object.
    * \param texImage is the target texture image.  It will have the texture \p
    * width, \p height, \p depth, \p border and \p internalFormat information.
    *      
    * \a retainInternalCopy is returned by this function and indicates whether
    * core Mesa should keep an internal copy of the texture image.
    */
   void (*CompressedTexImage1D)( MGLcontext *ctx, GLenum target,
                                 GLint level, GLint internalFormat,
                                 GLsizei width, GLint border,
                                 GLsizei imageSize, const GLvoid *data,
                                 struct gl_texture_object *texObj,
                                 struct gl_texture_image *texImage );
   /**
    * Called by glCompressedTexImage2D().
    *
    * \sa dd_function_table::CompressedTexImage1D.
    */
   void (*CompressedTexImage2D)( MGLcontext *ctx, GLenum target,
                                 GLint level, GLint internalFormat,
                                 GLsizei width, GLsizei height, GLint border,
                                 GLsizei imageSize, const GLvoid *data,
                                 struct gl_texture_object *texObj,
                                 struct gl_texture_image *texImage );
   /**
    * Called by glCompressedTexImage3D().
    *
    * \sa dd_function_table::CompressedTexImage3D.
    */
   void (*CompressedTexImage3D)( MGLcontext *ctx, GLenum target,
                                 GLint level, GLint internalFormat,
                                 GLsizei width, GLsizei height, GLsizei depth,
                                 GLint border,
                                 GLsizei imageSize, const GLvoid *data,
                                 struct gl_texture_object *texObj,
                                 struct gl_texture_image *texImage );

   /**
    * Called by glCompressedTexSubImage1D().
    * 
    * \param target user specified.
    * \param level user specified.
    * \param xoffset user specified.
    * \param yoffset user specified.
    * \param zoffset user specified.
    * \param width user specified.
    * \param height user specified.
    * \param depth user specified.
    * \param imageSize user specified.
    * \param data user specified.
    * \param texObj is the target texture object.
    * \param texImage is the target texture image.  It will have the texture \p
    * width, \p height, \p depth, \p border and \p internalFormat information.
    */
   void (*CompressedTexSubImage1D)(MGLcontext *ctx, GLenum target, GLint level,
                                   GLint xoffset, GLsizei width,
                                   GLenum format,
                                   GLsizei imageSize, const GLvoid *data,
                                   struct gl_texture_object *texObj,
                                   struct gl_texture_image *texImage);
   /**
    * Called by glCompressedTexSubImage2D().
    *
    * \sa dd_function_table::CompressedTexImage3D.
    */
   void (*CompressedTexSubImage2D)(MGLcontext *ctx, GLenum target, GLint level,
                                   GLint xoffset, GLint yoffset,
                                   GLsizei width, GLint height,
                                   GLenum format,
                                   GLsizei imageSize, const GLvoid *data,
                                   struct gl_texture_object *texObj,
                                   struct gl_texture_image *texImage);
   /**
    * Called by glCompressedTexSubImage3D().
    *
    * \sa dd_function_table::CompressedTexImage3D.
    */
   void (*CompressedTexSubImage3D)(MGLcontext *ctx, GLenum target, GLint level,
                                   GLint xoffset, GLint yoffset, GLint zoffset,
                                   GLsizei width, GLint height, GLint depth,
                                   GLenum format,
                                   GLsizei imageSize, const GLvoid *data,
                                   struct gl_texture_object *texObj,
                                   struct gl_texture_image *texImage);

   /**
    * Called to query number of bytes of storage needed to store the
    * specified compressed texture.
    */
   GLuint (*CompressedTextureSize)( MGLcontext *ctx, GLsizei width,
                                    GLsizei height, GLsizei depth,
                                    GLenum format );
   /*@}*/

   /**
    * \name Texture object functions
    */
   /*@{*/

   /**
    * Called by glBindTexture().
    */
   void (*BindTexture)( MGLcontext *ctx, GLenum target,
                        struct gl_texture_object *tObj );

   /**
    * Called to allocate a new texture object.
    * A new gl_texture_object should be returned.  The driver should
    * attach to it any device-specific info it needs.
    */
   struct gl_texture_object * (*NewTextureObject)( MGLcontext *ctx, GLuint name,
                                                   GLenum target );
   /**
    * Called when a texture object is about to be deallocated.  
    *
    * Driver should delete the gl_texture_object object and anything
    * hanging off of it.
    */
   void (*DeleteTexture)( MGLcontext *ctx, struct gl_texture_object *tObj );

   /**
    * Called to allocate a new texture image object.
    */
   struct gl_texture_image * (*NewTextureImage)( MGLcontext *ctx );

   /**
    * Called by glAreTextureResident().
    */
   GLboolean (*IsTextureResident)( MGLcontext *ctx,
                                   struct gl_texture_object *t );

   /**
    * Called by glPrioritizeTextures().
    */
   void (*PrioritizeTexture)( MGLcontext *ctx,  struct gl_texture_object *t,
                              GLclampf priority );

   /**
    * Called by glActiveTextureARB() to set current texture unit.
    */
   void (*ActiveTexture)( MGLcontext *ctx, GLuint texUnitNumber );

   /**
    * Called when the texture's color lookup table is changed.
    * 
    * If \p tObj is NULL then the shared texture palette
    * gl_texture_object::Palette is to be updated.
    */
   void (*UpdateTexturePalette)( MGLcontext *ctx,
                                 struct gl_texture_object *tObj );
   /*@}*/

   
   /**
    * \name Imaging functionality
    */
   /*@{*/
   void (*CopyColorTable)( MGLcontext *ctx,
			   GLenum target, GLenum internalformat,
			   GLint x, GLint y, GLsizei width );

   void (*CopyColorSubTable)( MGLcontext *ctx,
			      GLenum target, GLsizei start,
			      GLint x, GLint y, GLsizei width );

   void (*CopyConvolutionFilter1D)( MGLcontext *ctx, GLenum target,
				    GLenum internalFormat,
				    GLint x, GLint y, GLsizei width );

   void (*CopyConvolutionFilter2D)( MGLcontext *ctx, GLenum target,
				    GLenum internalFormat,
				    GLint x, GLint y,
				    GLsizei width, GLsizei height );
   /*@}*/


   /**
    * \name Vertex/fragment program functions
    */
   /*@{*/
   /** Bind a vertex/fragment program */
   void (*BindProgram)(MGLcontext *ctx, GLenum target, struct program *prog);
   /** Allocate a new program */
   struct program * (*NewProgram)(MGLcontext *ctx, GLenum target, GLuint id);
   /** Delete a program */
   void (*DeleteProgram)(MGLcontext *ctx, struct program *prog);   
   /** Notify driver that a program string has been specified. */
   void (*ProgramStringNotify)(MGLcontext *ctx, GLenum target, 
			       struct program *prog);
   


   /** Query if program can be loaded onto hardware */
   GLboolean (*IsProgramNative)(MGLcontext *ctx, GLenum target, 
				struct program *prog);
   
   /*@}*/


   /**
    * \name State-changing functions.
    *
    * \note drawing functions are above.
    *
    * These functions are called by their corresponding OpenGL API functions.
    * They are \e also called by the gl_PopAttrib() function!!!
    * May add more functions like these to the device driver in the future.
    */
   /*@{*/
   /** Specify the alpha test function */
   void (*AlphaFunc)(MGLcontext *ctx, GLenum func, GLfloat ref);
   /** Set the blend color */
   void (*BlendColor)(MGLcontext *ctx, const GLfloat color[4]);
   /** Set the blend equation */
   void (*BlendEquationSeparate)(MGLcontext *ctx, GLenum modeRGB, GLenum modeA);
   /** Specify pixel arithmetic */
   void (*BlendFuncSeparate)(MGLcontext *ctx,
                             GLenum sfactorRGB, GLenum dfactorRGB,
                             GLenum sfactorA, GLenum dfactorA);
   /** Specify clear values for the color buffers */
   void (*ClearColor)(MGLcontext *ctx, const GLfloat color[4]);
   /** Specify the clear value for the depth buffer */
   void (*ClearDepth)(MGLcontext *ctx, GLclampd d);
   /** Specify the clear value for the color index buffers */
   void (*ClearIndex)(MGLcontext *ctx, GLuint index);
   /** Specify the clear value for the stencil buffer */
   void (*ClearStencil)(MGLcontext *ctx, GLint s);
   /** Specify a plane against which all geometry is clipped */
   void (*ClipPlane)(MGLcontext *ctx, GLenum plane, const GLfloat *equation );
   /** Enable and disable writing of frame buffer color components */
   void (*ColorMask)(MGLcontext *ctx, GLboolean rmask, GLboolean gmask,
                     GLboolean bmask, GLboolean amask );
   /** Cause a material color to track the current color */
   void (*ColorMaterial)(MGLcontext *ctx, GLenum face, GLenum mode);
   /** Specify whether front- or back-facing facets can be culled */
   void (*CullFace)(MGLcontext *ctx, GLenum mode);
   /** Define front- and back-facing polygons */
   void (*FrontFace)(MGLcontext *ctx, GLenum mode);
   /** Specify the value used for depth buffer comparisons */
   void (*DepthFunc)(MGLcontext *ctx, GLenum func);
   /** Enable or disable writing into the depth buffer */
   void (*DepthMask)(MGLcontext *ctx, GLboolean flag);
   /** Specify mapping of depth values from NDC to window coordinates */
   void (*DepthRange)(MGLcontext *ctx, GLclampd nearval, GLclampd farval);
   /** Specify the current buffer for writing */
   void (*DrawBuffer)( MGLcontext *ctx, GLenum buffer );
   /** Enable or disable server-side gl capabilities */
   void (*Enable)(MGLcontext *ctx, GLenum cap, GLboolean state);
   /** Specify fog parameters */
   void (*Fogfv)(MGLcontext *ctx, GLenum pname, const GLfloat *params);
   /** Specify implementation-specific hints */
   void (*Hint)(MGLcontext *ctx, GLenum target, GLenum mode);
   /** Control the writing of individual bits in the color index buffers */
   void (*IndexMask)(MGLcontext *ctx, GLuint mask);
   /** Set light source parameters */
   void (*Lightfv)(MGLcontext *ctx, GLenum light,
		   GLenum pname, const GLfloat *params );
   /** Set the lighting model parameters */
   void (*LightModelfv)(MGLcontext *ctx, GLenum pname, const GLfloat *params);
   /** Specify the line stipple pattern */
   void (*LineStipple)(MGLcontext *ctx, GLint factor, GLushort pattern );
   /** Specify the width of rasterized lines */
   void (*LineWidth)(MGLcontext *ctx, GLfloat width);
   /** Specify a logical pixel operation for color index rendering */
   void (*LogicOpcode)(MGLcontext *ctx, GLenum opcode);
   void (*PointParameterfv)(MGLcontext *ctx, GLenum pname,
                            const GLfloat *params);
   /** Specify the diameter of rasterized points */
   void (*PointSize)(MGLcontext *ctx, GLfloat size);
   /** Select a polygon rasterization mode */
   void (*PolygonMode)(MGLcontext *ctx, GLenum face, GLenum mode);
   /** Set the scale and units used to calculate depth values */
   void (*PolygonOffset)(MGLcontext *ctx, GLfloat factor, GLfloat units);
   /** Set the polygon stippling pattern */
   void (*PolygonStipple)(MGLcontext *ctx, const GLubyte *mask );
   /* Specifies the current buffer for reading */
   void (*ReadBuffer)( MGLcontext *ctx, GLenum buffer );
   /** Set rasterization mode */
   void (*RenderMode)(MGLcontext *ctx, GLenum mode );
   /** Define the scissor box */
   void (*Scissor)(MGLcontext *ctx, GLint x, GLint y, GLsizei w, GLsizei h);
   /** Select flat or smooth shading */
   void (*ShadeModel)(MGLcontext *ctx, GLenum mode);
   /** Set function and reference value for stencil testing */
   void (*StencilFunc)(MGLcontext *ctx, GLenum func, GLint ref, GLuint mask);
   /** Control the writing of individual bits in the stencil planes */
   void (*StencilMask)(MGLcontext *ctx, GLuint mask);
   /** Set stencil test actions */
   void (*StencilOp)(MGLcontext *ctx, GLenum fail, GLenum zfail, GLenum zpass);
   void (*ActiveStencilFace)(MGLcontext *ctx, GLuint face);
   /** Control the generation of texture coordinates */
   void (*TexGen)(MGLcontext *ctx, GLenum coord, GLenum pname,
		  const GLfloat *params);
   /** Set texture environment parameters */
   void (*TexEnv)(MGLcontext *ctx, GLenum target, GLenum pname,
                  const GLfloat *param);
   /** Set texture parameters */
   void (*TexParameter)(MGLcontext *ctx, GLenum target,
                        struct gl_texture_object *texObj,
                        GLenum pname, const GLfloat *params);
   void (*TextureMatrix)(MGLcontext *ctx, GLuint unit, const MGLmatrix *mat);
   /** Set the viewport */
   void (*Viewport)(MGLcontext *ctx, GLint x, GLint y, GLsizei w, GLsizei h);
   /*@}*/


   /**
    * \name Vertex array functions
    *
    * Called by the corresponding OpenGL functions.
    */
   /*@{*/
   void (*VertexPointer)(MGLcontext *ctx, GLint size, GLenum type,
			 GLsizei stride, const GLvoid *ptr);
   void (*NormalPointer)(MGLcontext *ctx, GLenum type,
			 GLsizei stride, const GLvoid *ptr);
   void (*ColorPointer)(MGLcontext *ctx, GLint size, GLenum type,
			GLsizei stride, const GLvoid *ptr);
   void (*FogCoordPointer)(MGLcontext *ctx, GLenum type,
			   GLsizei stride, const GLvoid *ptr);
   void (*IndexPointer)(MGLcontext *ctx, GLenum type,
			GLsizei stride, const GLvoid *ptr);
   void (*SecondaryColorPointer)(MGLcontext *ctx, GLint size, GLenum type,
				 GLsizei stride, const GLvoid *ptr);
   void (*TexCoordPointer)(MGLcontext *ctx, GLint size, GLenum type,
			   GLsizei stride, const GLvoid *ptr);
   void (*EdgeFlagPointer)(MGLcontext *ctx, GLsizei stride, const GLvoid *ptr);
   void (*VertexAttribPointer)(MGLcontext *ctx, GLuint index, GLint size,
                               GLenum type, GLsizei stride, const GLvoid *ptr);
   void (*LockArraysEXT)( MGLcontext *ctx, GLint first, GLsizei count );
   void (*UnlockArraysEXT)( MGLcontext *ctx );
   /*@}*/


   /** 
    * \name State-query functions
    *
    * Return GL_TRUE if query was completed, GL_FALSE otherwise.
    */
   /*@{*/
   /** Return the value or values of a selected parameter */
   GLboolean (*GetBooleanv)(MGLcontext *ctx, GLenum pname, GLboolean *result);
   /** Return the value or values of a selected parameter */
   GLboolean (*GetDoublev)(MGLcontext *ctx, GLenum pname, GLdouble *result);
   /** Return the value or values of a selected parameter */
   GLboolean (*GetFloatv)(MGLcontext *ctx, GLenum pname, GLfloat *result);
   /** Return the value or values of a selected parameter */
   GLboolean (*GetIntegerv)(MGLcontext *ctx, GLenum pname, GLint *result);
   /** Return the value or values of a selected parameter */
   GLboolean (*GetPointerv)(MGLcontext *ctx, GLenum pname, GLvoid **result);
   /*@}*/
   

   /**
    * \name Vertex buffer object functions
    */
#if FEATURE_ARB_vertex_buffer_object
   /*@{*/
   void (*BindBuffer)( MGLcontext *ctx, GLenum target,
		       struct Mgl_buffer_object *obj );

   struct Mgl_buffer_object * (*NewBufferObject)( MGLcontext *ctx, GLuint buffer,
						 GLenum target );
   
   void (*DeleteBuffer)( MGLcontext *ctx, struct Mgl_buffer_object *obj );

   void (*BufferData)( MGLcontext *ctx, GLenum target, GLsizeiptrARB size,
		       const GLvoid *data, GLenum usage,
		       struct Mgl_buffer_object *obj );

   void (*BufferSubData)( MGLcontext *ctx, GLenum target, GLintptrARB offset,
			  GLsizeiptrARB size, const GLvoid *data,
			  struct Mgl_buffer_object *obj );

   void (*GetBufferSubData)( MGLcontext *ctx, GLenum target,
			     GLintptrARB offset, GLsizeiptrARB size,
			     GLvoid *data, struct Mgl_buffer_object *obj );

   void * (*MapBuffer)( MGLcontext *ctx, GLenum target, GLenum access,
			struct Mgl_buffer_object *obj );

   GLboolean (*UnmapBuffer)( MGLcontext *ctx, GLenum target,
			     struct Mgl_buffer_object *obj );
   /*@}*/
#endif

   /**
    * \name Support for multiple T&L engines
    */
   /*@{*/

   /**
    * Bitmask of state changes that require the current T&L module to be
    * validated, using ValidateTnlModule() below.
    */
   GLuint NeedValidate;

   /**
    * Validate the current T&L module. 
    *
    * This is called directly after UpdateState() when a state change that has
    * occurred matches the dd_function_table::NeedValidate bitmask above.  This
    * ensures all computed values are up to date, thus allowing the driver to
    * decide if the current T&L module needs to be swapped out.
    *
    * This must be non-NULL if a driver installs a custom T&L module and sets
    * the dd_function_table::NeedValidate bitmask, but may be NULL otherwise.
    */
   void (*ValidateTnlModule)( MGLcontext *ctx, GLuint new_state );


#define PRIM_OUTSIDE_BEGIN_END   GL_POLYGON+1
#define PRIM_INSIDE_UNKNOWN_PRIM GL_POLYGON+2
#define PRIM_UNKNOWN             GL_POLYGON+3

   /**
    * Set by the driver-supplied T&L engine.  
    *
    * Set to PRIM_OUTSIDE_BEGIN_END when outside glBegin()/glEnd().
    */
   GLuint CurrentExecPrimitive;

   /**
    * Current state of an in-progress compilation.  
    *
    * May take on any of the additional values PRIM_OUTSIDE_BEGIN_END,
    * PRIM_INSIDE_UNKNOWN_PRIM or PRIM_UNKNOWN defined above.
    */
   GLuint CurrentSavePrimitive;


#define FLUSH_STORED_VERTICES 0x1
#define FLUSH_UPDATE_CURRENT  0x2
   /**
    * Set by the driver-supplied T&L engine whenever vertices are buffered
    * between glBegin()/glEnd() objects or __GLcontextRec::Current is not
    * updated.
    *
    * The dd_function_table::FlushVertices call below may be used to resolve
    * these conditions.
    */
   GLuint NeedFlush;
   GLuint SaveNeedFlush;

   /**
    * If inside glBegin()/glEnd(), it should ASSERT(0).  Otherwise, if
    * FLUSH_STORED_VERTICES bit in \p flags is set flushes any buffered
    * vertices, if FLUSH_UPDATE_CURRENT bit is set updates
    * __GLcontextRec::Current and gl_light_attrib::Material
    *
    * Note that the default T&L engine never clears the
    * FLUSH_UPDATE_CURRENT bit, even after performing the update.
    */
   void (*FlushVertices)( MGLcontext *ctx, GLuint flags );
   void (*SaveFlushVertices)( MGLcontext *ctx );

   /**
    * Give the driver the opportunity to hook in its own vtxfmt for
    * compiling optimized display lists.  This is called on each valid
    * glBegin() during list compilation.
    */
   GLboolean (*NotifySaveBegin)( MGLcontext *ctx, GLenum mode );

   /**
    * Notify driver that the special derived value _NeedEyeCoords has
    * changed.
    */
   void (*LightingSpaceChange)( MGLcontext *ctx );

   /**
    * Let the T&L component know when the context becomes current.
    */
   void (*MakeCurrent)( MGLcontext *ctx, MGLframebuffer *drawBuffer, MGLframebuffer *readBuffer );

   /**
    * Called by glNewList().
    *
    * Let the T&L component know what is going on with display lists
    * in time to make changes to dispatch tables, etc.
    */
   void (*NewList)( MGLcontext *ctx, GLuint list, GLenum mode );
   /**
    * Called by glEndList().
    *
    * \sa dd_function_table::NewList.
    */
   void (*EndList)( MGLcontext *ctx );

   /**
    * Called by glCallList(s), but not recursively.
    *
    * Notify the T&L component before and after calling a display list.
    * Called by glCallList(s), but not recursively.
    */
   void (*BeginCallList)( MGLcontext *ctx, GLuint list );
   /**
    * Called by glEndCallList().
    *
    * \sa dd_function_table::BeginCallList.
    */
   void (*EndCallList)( MGLcontext *ctx );
};

/**
 * Constants which may be overridden by device driver during context creation
 * but are never changed after that.
 */
struct Mgl_constants
{
   GLint MaxTextureLevels;		/**< Maximum number of allowed mipmap levels. */ 
   GLint Max3DTextureLevels;		/**< Maximum number of allowed mipmap levels for 3D texture targets. */
   GLint MaxCubeTextureLevels;          /**< Maximum number of allowed mipmap levels for GL_ARB_texture_cube_map */
   GLint MaxTextureRectSize;            /* GL_NV_texture_rectangle */
   GLuint MaxTextureCoordUnits;
   GLuint MaxTextureImageUnits;
   GLuint MaxTextureUnits;              /* = MAX(CoordUnits, ImageUnits) */
   GLfloat MaxTextureMaxAnisotropy;	/* GL_EXT_texture_filter_anisotropic */
   GLfloat MaxTextureLodBias;           /* GL_EXT_texture_lod_bias */
   GLuint MaxArrayLockSize;
   GLint SubPixelBits;
   GLfloat MinPointSize, MaxPointSize;		/* aliased */
   GLfloat MinPointSizeAA, MaxPointSizeAA;	/* antialiased */
   GLfloat PointSizeGranularity;
   GLfloat MinLineWidth, MaxLineWidth;		/* aliased */
   GLfloat MinLineWidthAA, MaxLineWidthAA;	/* antialiased */
   GLfloat LineWidthGranularity;
   GLuint MaxColorTableSize;
   GLuint MaxConvolutionWidth;
   GLuint MaxConvolutionHeight;
   GLuint MaxClipPlanes;
   GLuint MaxLights;
   GLfloat MaxShininess;			/* GL_NV_light_max_exponent */
   GLfloat MaxSpotExponent;			/* GL_NV_light_max_exponent */
   /* GL_ARB_vertex_program */
   GLuint MaxVertexProgramInstructions;
   GLuint MaxVertexProgramAttribs;
   GLuint MaxVertexProgramTemps;
   GLuint MaxVertexProgramLocalParams;
   GLuint MaxVertexProgramEnvParams;
   GLuint MaxVertexProgramAddressRegs;
   /* GL_ARB_fragment_program */
   GLuint MaxFragmentProgramInstructions;
   GLuint MaxFragmentProgramAttribs;
   GLuint MaxFragmentProgramTemps;
   GLuint MaxFragmentProgramLocalParams;
   GLuint MaxFragmentProgramEnvParams;
   GLuint MaxFragmentProgramAddressRegs;
   GLuint MaxFragmentProgramAluInstructions;
   GLuint MaxFragmentProgramTexInstructions;
   GLuint MaxFragmentProgramTexIndirections;
   /* vertex or fragment program */
   GLuint MaxProgramMatrices;
   GLuint MaxProgramMatrixStackDepth;
   /* vertex array / buffer object bounds checking */
   GLboolean CheckArrayBounds;
};

/**
 * Mesa context
 *
 * This is the central context data structure for Mesa.  Almost all
 * OpenGL state is contained in this structure.
 * Think of this as a base class from which device drivers will derive
 * sub classes.
 */
struct M__GLcontextRec {
   /**
    * \name OS related interfaces. 
    *
    * These \b must be the first members of this structure, because they are
    * exposed to the outside world (i.e. GLX extension).
    */
   /*@{*/
   //__GLimports imports; //OC
   //__GLexports exports;
   /*@}*/

   /** State possibly shared with other contexts in the address space */
   struct gl_shared_state *Shared;

   /** \name API function pointer tables */
   /*@{*/
   struct _glapi_table *Save;	/**< Display list save functions */
   struct _glapi_table *Exec;	/**< Execute functions */
   struct _glapi_table *CurrentDispatch;  /**< == Save or Exec !! */
   /*@}*/

   MGLvisual Visual;
   MGLframebuffer *DrawBuffer;	/**< buffer for writing */
   MGLframebuffer *ReadBuffer;	/**< buffer for reading */

   /**
    * Device driver function pointer table
    */
   struct Mdd_function_table Driver;

   void *DriverCtx;	/**< Points to device driver context/state */
   void *DriverMgrCtx;	/**< Points to device driver manager (optional)*/

   /** Core/Driver constants */
   struct Mgl_constants Const;

   /** \name The various 4x4 matrix stacks */
   /*@{*/
   struct Mmatrix_stack ModelviewMatrixStack;
   struct Mmatrix_stack ProjectionMatrixStack;
   struct Mmatrix_stack ColorMatrixStack;
   struct Mmatrix_stack TextureMatrixStack[MMAX_TEXTURE_COORD_UNITS];
   struct Mmatrix_stack ProgramMatrixStack[MMAX_PROGRAM_MATRICES];
   struct Mmatrix_stack *CurrentStack; /**< Points to one of the above stacks */
   /*@}*/

   /** Combined modelview and projection matrix */
   MGLmatrix _ModelProjectMatrix;

   /** \name Display lists */
   struct Mmesa_list_state ListState;

   GLboolean ExecuteFlag;	/**< Execute GL commands? */
   GLboolean CompileFlag;	/**< Compile GL commands into display list? */

   /** Extensions */
   struct Mgl_extensions Extensions;

   /** \name Renderer attribute stack */
   /*@{*/
   GLuint AttribStackDepth;
   struct gl_attrib_node *AttribStack[MMAX_ATTRIB_STACK_DEPTH];
   /*@}*/

   /** \name Renderer attribute groups
    * 
    * We define a struct for each attribute group to make pushing and popping
    * attributes easy.  Also it's a good organization.
    */
   /*@{*/
   struct Mgl_accum_attrib	Accum;		/**< Accumulation buffer attributes */
   struct Mgl_colorbuffer_attrib	Color;		/**< Color buffers attributes */
   struct Mgl_current_attrib	Current;	/**< Current attributes */
   struct Mgl_depthbuffer_attrib	Depth;		/**< Depth buffer attributes */
   struct Mgl_eval_attrib	Eval;		/**< Eval attributes */
   struct Mgl_fog_attrib		Fog;		/**< Fog attributes */
   struct Mgl_hint_attrib	Hint;		/**< Hint attributes */
   struct Mgl_light_attrib	Light;		/**< Light attributes */
   struct Mgl_line_attrib	Line;		/**< Line attributes */
   struct Mgl_list_attrib	List;		/**< List attributes */
   struct Mgl_multisample_attrib Multisample;
   struct Mgl_pixel_attrib	Pixel;		/**< Pixel attributes */
   struct Mgl_point_attrib	Point;		/**< Point attributes */
   struct Mgl_polygon_attrib	Polygon;	/**< Polygon attributes */
   GLuint PolygonStipple[32];			/**< Polygon stipple */
   struct Mgl_scissor_attrib	Scissor;	/**< Scissor attributes */
   struct Mgl_stencil_attrib	Stencil;	/**< Stencil buffer attributes */
   struct Mgl_texture_attrib	Texture;	/**< Texture attributes */
   struct Mgl_transform_attrib	Transform;	/**< Transformation attributes */
   struct Mgl_viewport_attrib	Viewport;	/**< Viewport attributes */
   /*@}*/

   /** \name Other attribute groups */
   /*@{*/
   struct Mgl_histogram_attrib	Histogram;
   struct Mgl_minmax_attrib	MinMax;
   struct Mgl_convolution_attrib Convolution1D;
   struct Mgl_convolution_attrib Convolution2D;
   struct Mgl_convolution_attrib Separable2D;
   /*@}*/

   /** \name Client attribute stack */
   /*@{*/
   GLuint ClientAttribStackDepth;
   struct gl_attrib_node *ClientAttribStack[MMAX_CLIENT_ATTRIB_STACK_DEPTH];
   /*@}*/

   /** \name Client attribute groups */
   /*@{*/
   struct Mgl_array_attrib	Array;	/**< Vertex arrays */
   struct Mgl_pixelstore_attrib	Pack;	/**< Pixel packing */
   struct Mgl_pixelstore_attrib	Unpack;	/**< Pixel unpacking */
   struct Mgl_pixelstore_attrib	DefaultPacking;	/**< Default params */

   struct Mgl_evaluators EvalMap;   /**< All evaluators */
   struct Mgl_feedback   Feedback;  /**< Feedback */
   struct Mgl_selection  Select;    /**< Selection */

   struct Mgl_color_table ColorTable;       /**< Pre-convolution */
   struct Mgl_color_table ProxyColorTable;  /**< Pre-convolution */
   struct Mgl_color_table PostConvolutionColorTable;
   struct Mgl_color_table ProxyPostConvolutionColorTable;
   struct Mgl_color_table PostColorMatrixColorTable;
   struct Mgl_color_table ProxyPostColorMatrixColorTable;

   struct Mprogram_state Program;             /**< for vertex or fragment progs */
   struct Mvertex_program_state VertexProgram;      /**< GL_NV_vertex_program */
   struct Mfragment_program_state FragmentProgram;  /**< GL_NV_fragment_program */

   struct Mocclusion_state Occlusion;  /**< GL_ARB_occlusion_query */

   GLenum ErrorValue;        /**< Last error code */
   GLenum RenderMode;        /**< either GL_RENDER, GL_SELECT, GL_FEEDBACK */
   GLuint NewState;          /**< bitwise-or of _NEW_* flags */
   /*@}*/

   /** \name Derived */
   /*@{*/
   GLuint _TriangleCaps;      /**< bitwise-or of DD_* flags */
   GLuint _ImageTransferState;/**< bitwise-or of IMAGE_*_BIT flags */
   GLfloat _EyeZDir[3];
   GLfloat _ModelViewInvScale;
   GLuint _NeedEyeCoords;
   GLuint _ForceEyeCoords; 
   GLboolean _RotateMode;
   GLenum _CurrentProgram;    /* currently executing program */

   struct Mgl_shine_tab *_ShineTable[2]; /**< Active shine tables */
   struct Mgl_shine_tab *_ShineTabList;  /**< MRU list of inactive shine tables */
   /**@}*/

   struct Mgl_list_extensions ListExt; /**< driver dlist extensions */


   GLboolean OcclusionResult;       /**< for GL_HP_occlusion_test */
   GLboolean OcclusionResultSaved;  /**< for GL_HP_occlusion_test */
   GLuint _Facing; /**< This is a hack for 2-sided stencil test.
		    *
		    * We don't have a better way to communicate this value from
		    * swrast_setup to swrast. */


   /** \name Z buffer stuff */
   /*@{*/
   GLuint DepthMax;	/**< Max depth buffer value */
   GLfloat DepthMaxF;	/**< Float max depth buffer value */
   GLfloat MRD;		/**< minimum resolvable difference in Z values */
   /*@}*/

   /** \name Color clamping (tentative part of GL_ARB_color_clamp_control) */
   /*@{*/
   GLboolean ClampFragmentColors;
   GLboolean ClampVertexColors;
   /*@}*/

   /** \name For debugging/development only */
   /*@{*/
   GLboolean FirstTimeCurrent;
   /*@}*/

   /** Dither disable via MESA_NO_DITHER env var */
   GLboolean NoDither;

   /** Core tnl module support */
   struct Mgl_tnl_module TnlModule;

   /**
    * \name Hooks for module contexts.  
    *
    * These will eventually live in the driver or elsewhere.
    */
   /*@{*/
   void *swrast_context;
   void *swsetup_context;
   void *swtnl_context;
   void *swtnl_im;
   void *acache_context;
   void *aelt_context;
   /*@}*/
};

struct Mosmesa_context {
   MGLcontext mesa;		/* The core GL/Mesa context */
   MGLvisual *gl_visual;		/* Describes the buffers */
   MGLframebuffer *gl_buffer;	/* Depth, stencil, accum, etc buffers */
   GLenum format;		/* either GL_RGBA or GL_COLOR_INDEX */
   void *buffer;		/* the image buffer */
   GLint width, height;		/* size of image buffer */
   GLint rowlength;		/* number of pixels per row */
   GLint userRowLength;		/* user-specified number of pixels per row */
   GLint rshift, gshift;	/* bit shifts for RGBA formats */
   GLint bshift, ashift;
   GLint rInd, gInd, bInd, aInd;/* index offsets for RGBA formats */
   MGLchan *rowaddr[MMAX_HEIGHT];	/* address of first pixel in each image row */
   GLboolean yup;		/* TRUE  -> Y increases upward */
				/* FALSE -> Y increases downward */
};

typedef struct {
   GLuint NewState;
   GLenum render_prim;
   GLuint last_index;
   //SWvertex *verts;
} MSScontext;

#define MSWSETUP_CONTEXT(ctx) ((MSScontext *)ctx->swsetup_context)

//-------------------------------------------------------------------------

/*
 * Note: The first attributes match the VERT_ATTRIB_* definitions
 * in mtypes.h.  However, the tnl module has additional attributes
 * for materials, color indexes, edge flags, etc.
 */
/* Note: These are currently being used to define both inputs and
 * outputs from the tnl pipeline.  A better solution (which would also
 * releive the congestion to slightly prolong the life of the bitmask
 * below) is to have the fixed function pipeline populate a set of
 * arrays named after those produced by the vertex program stage, and
 * have the rest the mesa backend work on those.
 */
enum {
	M_TNL_ATTRIB_POS = 0,
	M_TNL_ATTRIB_WEIGHT = 1,
	M_TNL_ATTRIB_NORMAL = 2,
	M_TNL_ATTRIB_COLOR0 = 3,
	M_TNL_ATTRIB_COLOR1 = 4,
	M_TNL_ATTRIB_FOG = 5,
	M_TNL_ATTRIB_SIX = 6,
	M_TNL_ATTRIB_SEVEN = 7,
	M_TNL_ATTRIB_TEX0 = 8,
	M_TNL_ATTRIB_TEX1 = 9,
	M_TNL_ATTRIB_TEX2 = 10,
	M_TNL_ATTRIB_TEX3 = 11,
	M_TNL_ATTRIB_TEX4 = 12,
	M_TNL_ATTRIB_TEX5 = 13,
	M_TNL_ATTRIB_TEX6 = 14,
	M_TNL_ATTRIB_TEX7 = 15,
	M_TNL_ATTRIB_MAT_FRONT_AMBIENT = 16,
	M_TNL_ATTRIB_MAT_BACK_AMBIENT = 17,
	M_TNL_ATTRIB_MAT_FRONT_DIFFUSE = 18,
	M_TNL_ATTRIB_MAT_BACK_DIFFUSE = 19,
	M_TNL_ATTRIB_MAT_FRONT_SPECULAR = 20,
	M_TNL_ATTRIB_MAT_BACK_SPECULAR = 21,
	M_TNL_ATTRIB_MAT_FRONT_EMISSION = 22,
	M_TNL_ATTRIB_MAT_BACK_EMISSION = 23,
	M_TNL_ATTRIB_MAT_FRONT_SHININESS = 24,
	M_TNL_ATTRIB_MAT_BACK_SHININESS = 25,
	M_TNL_ATTRIB_MAT_FRONT_INDEXES = 26,
	M_TNL_ATTRIB_MAT_BACK_INDEXES = 27, 
	M_TNL_ATTRIB_INDEX = 28,        
	M_TNL_ATTRIB_EDGEFLAG = 29,     
	M_TNL_ATTRIB_POINTSIZE = 30,
	M_TNL_ATTRIB_MAX = 31
} ;

struct Mtnl_clipspace;
struct Mtnl_clipspace_attr;

typedef void (*Mtnl_extract_func)( const struct Mtnl_clipspace_attr *a, GLfloat *out, const GLubyte *v );
typedef void (*Mtnl_insert_func)( const struct Mtnl_clipspace_attr *a, GLubyte *v, const GLfloat *in );
typedef void (*Mtnl_emit_func)( MGLcontext *ctx, GLuint start, GLuint end, void *dest );

typedef void (*Mtnl_points_func)( MGLcontext *ctx, GLuint first, GLuint last );
typedef void (*Mtnl_line_func)( MGLcontext *ctx, GLuint v1, GLuint v2 );
typedef void (*Mtnl_triangle_func)( MGLcontext *ctx, GLuint v1, GLuint v2, GLuint v3 );
typedef void (*Mtnl_quad_func)( MGLcontext *ctx, GLuint v1, GLuint v2, GLuint v3, GLuint v4 );
typedef void (*Mtnl_render_func)( MGLcontext *ctx, GLuint start, GLuint count, GLuint flags );
typedef void (*Mtnl_interp_func)( MGLcontext *ctx, GLfloat t, GLuint dst, GLuint out, GLuint in, GLboolean force_boundary );
typedef void (*Mtnl_copy_pv_func)( MGLcontext *ctx, GLuint dst, GLuint src );
typedef void (*Mtnl_setup_func)( MGLcontext *ctx, GLuint start, GLuint end, GLuint new_inputs);

typedef void (*Mtnl_attrfv_func)( const GLfloat * );

struct M_tnl_dynfn {
   struct M_tnl_dynfn *next, *prev;
   GLuint key;
   char *code;
};

struct Mtnl_clipspace_codegen {
   GLboolean (*emit_header)( struct tnl_clipspace_codegen *, struct tnl_clipspace *);
   GLboolean (*emit_footer)( struct tnl_clipspace_codegen * );
   GLboolean (*emit_attr_header)( struct tnl_clipspace_codegen *, struct tnl_clipspace_attr *, GLint j, GLenum out_type, GLboolean need_vp );
   GLboolean (*emit_attr_footer)( struct tnl_clipspace_codegen * );
   GLboolean (*emit_mov)( struct tnl_clipspace_codegen *, GLint, GLint );
   GLboolean (*emit_const)( struct tnl_clipspace_codegen *, GLint, GLfloat );
   GLboolean (*emit_mad)( struct tnl_clipspace_codegen *, GLint, GLint, GLint, GLint );
   GLboolean (*emit_float_to_chan)( struct tnl_clipspace_codegen *, GLint, GLint );
   GLboolean (*emit_const_chan)( struct tnl_clipspace_codegen *, GLint, MGLchan );
   GLboolean (*emit_float_to_ubyte)( struct tnl_clipspace_codegen *, GLint, GLint );
   GLboolean (*emit_const_ubyte)( struct tnl_clipspace_codegen *, GLint, GLubyte );
   Mtnl_emit_func (*emit_store_func)( struct tnl_clipspace_codegen * );
   
   struct M_tnl_dynfn codegen_list;
   
   char *buf;
   int buf_size;
   int buf_used;
   int out_offset;
};

/**
 * Describes how to convert/move a vertex attribute from a vertex array
 * to a vertex structure.
 */
struct Mtnl_clipspace_attr
{
   GLuint attrib;          /* which vertex attrib (0=position, etc) */
   GLuint format;
   GLuint vertoffset;      /* position of the attrib in the vertex struct */
   GLuint vertattrsize;    /* size of the attribute in bytes */
   GLubyte *inputptr;
   GLuint inputstride;
   Mtnl_insert_func *insert;
   Mtnl_insert_func emit;
   Mtnl_extract_func extract;
   const GLfloat *vp;   /* NDC->Viewport mapping matrix */
};

/**
 * Used to describe conversion of vertex arrays to vertex structures.
 * I.e. Structure of arrays to arrays of structs.
 */
struct Mtnl_clipspace
{
   GLboolean need_extras;
   
   GLuint new_inputs;

   GLubyte *vertex_buf;
   GLuint vertex_size;
   GLuint max_vertex_size;

   struct Mtnl_clipspace_attr attr[M_TNL_ATTRIB_MAX];
   GLuint attr_count;

   Mtnl_emit_func emit;
   Mtnl_interp_func interp;
   Mtnl_copy_pv_func copy_pv;

   struct Mtnl_clipspace_codegen codegen;
};

struct Mtnl_device_driver
{
   /***
    *** TNL Pipeline
    ***/

   void (*RunPipeline)(MGLcontext *ctx);
   /* Replaces PipelineStart/PipelineFinish -- intended to allow
    * drivers to wrap _tnl_run_pipeline() with code to validate state
    * and grab/release hardware locks.  
    */

   void (*NotifyMaterialChange)(MGLcontext *ctx);
   /* Alert tnl-aware drivers of changes to material.
    */

   GLboolean (*NotifyBegin)(MGLcontext *ctx, GLenum p);
   /* Allow drivers to hook in optimized begin/end engines.
    * Return value:  GL_TRUE - driver handled the begin
    *                GL_FALSE - driver didn't handle the begin
    */

   /***
    *** Rendering -- These functions called only from t_vb_render.c
    ***/
   struct
   {
      void (*Start)(MGLcontext *ctx);
      void (*Finish)(MGLcontext *ctx);
      /* Called before and after all rendering operations, including DrawPixels,
       * ReadPixels, Bitmap, span functions, and CopyTexImage, etc commands.
       * These are a suitable place for grabbing/releasing hardware locks.
       */

      void (*PrimitiveNotify)(MGLcontext *ctx, GLenum mode);
      /* Called between RenderStart() and RenderFinish() to indicate the
       * type of primitive we're about to draw.  Mode will be one of the
       * modes accepted by glBegin().
       */

      Mtnl_interp_func Interp;
      /* The interp function is called by the clipping routines when we need
       * to generate an interpolated vertex.  All pertinant vertex ancilliary
       * data should be computed by interpolating between the 'in' and 'out'
       * vertices.
       */

      Mtnl_copy_pv_func CopyPV;
      /* The copy function is used to make a copy of a vertex.  All pertinant
       * vertex attributes should be copied.
       */

      void (*ClippedPolygon)( MGLcontext *ctx, const GLuint *elts, GLuint n );
      /* Render a polygon with <n> vertices whose indexes are in the <elts>
       * array.
       */

      void (*ClippedLine)( MGLcontext *ctx, GLuint v0, GLuint v1 );
      /* Render a line between the two vertices given by indexes v0 and v1. */

      Mtnl_points_func           Points; /* must now respect vb->elts */
      Mtnl_line_func             Line;
      Mtnl_triangle_func         Triangle;
      Mtnl_quad_func             Quad;
      /* These functions are called in order to render points, lines,
       * triangles and quads.  These are only called via the T&L module.
       */

      Mtnl_render_func          *PrimTabVerts;
      Mtnl_render_func          *PrimTabElts;
      /* Render whole unclipped primitives (points, lines, linestrips,
       * lineloops, etc).  The tables are indexed by the GL enum of the
       * primitive to be rendered.  RenderTabVerts is used for non-indexed
       * arrays of vertices.  RenderTabElts is used for indexed arrays of
       * vertices.
       */

      void (*ResetLineStipple)( MGLcontext *ctx );
      /* Reset the hardware's line stipple counter.
       */

      Mtnl_setup_func BuildVertices;
      /* This function is called whenever new vertices are required for
       * rendering.  The vertices in question are those n such that start
       * <= n < end.  The new_inputs parameter indicates those fields of
       * the vertex which need to be updated, if only a partial repair of
       * the vertex is required.
       *
       * This function is called only from _tnl_render_stage in tnl/t_render.c.
       */
      

      GLboolean (*Multipass)( MGLcontext *ctx, GLuint passno );
      /* Driver may request additional render passes by returning GL_TRUE
       * when this function is called.  This function will be called
       * after the first pass, and passes will be made until the function
       * returns GL_FALSE.  If no function is registered, only one pass
       * is made.
       *
       * This function will be first invoked with passno == 1.
       */
   } Render;
};

struct Mtnl_prim {
   GLuint mode;
   GLuint start;
   GLuint count;
};

struct Mtnl_copied_vtx {
   GLfloat buffer[M_TNL_ATTRIB_MAX * 4 * MTNL_MAX_COPIED_VERTS];
   GLuint nr;
};

struct M_tnl_dynfn_lists {
   struct M_tnl_dynfn Vertex[4];
   struct M_tnl_dynfn Attribute[4];
};

struct M_tnl_dynfn_generators {
   struct M_tnl_dynfn *(*Vertex[4])( MGLcontext *ctx, int key );
   struct M_tnl_dynfn *(*Attribute[4])( MGLcontext *ctx, int key );
};

struct Mtnl_eval1_map {
   struct Mgl_1d_map *map;
   GLuint sz;
};

struct Mtnl_eval2_map {
   struct Mgl_2d_map *map;
   GLuint sz;
};

struct Mtnl_eval {
   GLuint new_state;
   struct Mtnl_eval1_map map1[M_TNL_ATTRIB_INDEX + 1];
   struct Mtnl_eval2_map map2[M_TNL_ATTRIB_INDEX + 1];
};

/* The assembly of vertices in immediate mode is separated from
 * display list compilation.  This allows a simpler immediate mode
 * treatment and a display list compiler better suited to
 * hardware-acceleration.
 */
struct Mtnl_vtx {
   GLfloat buffer[MVERT_BUFFER_SIZE];
   GLubyte attrsz[M_TNL_ATTRIB_MAX];
   GLuint vertex_size;
   struct Mtnl_prim prim[MTNL_MAX_PRIM];
   GLuint prim_count;
   GLfloat *vbptr;		      /* cursor, points into buffer */
   GLfloat vertex[M_TNL_ATTRIB_MAX*4]; /* current vertex */
   GLfloat *attrptr[M_TNL_ATTRIB_MAX]; /* points into vertex */
   GLfloat *current[M_TNL_ATTRIB_MAX]; /* points into ctx->Current, etc */
   GLuint counter, initial_counter;
   struct Mtnl_copied_vtx copied;

   Mtnl_attrfv_func tabfv[M_TNL_MAX_ATTR_CODEGEN+1][4]; /* plus 1 for ERROR_ATTRIB */

   struct M_tnl_dynfn_lists cache;
   struct M_tnl_dynfn_generators gen;

   struct Mtnl_eval eval;
   GLboolean *edgeflag_tmp;
   GLboolean have_materials;
};

struct Mtnl_save {
   GLubyte attrsz[M_TNL_ATTRIB_MAX];
   GLuint vertex_size;

   GLfloat *buffer;
   GLuint count;
   GLuint wrap_count;

   struct tnl_prim *prim;
   GLuint prim_count, prim_max;

   struct Mtnl_vertex_store *vertex_store;
   struct Mtnl_primitive_store *prim_store;

   GLfloat *vbptr;		   /* cursor, points into buffer */
   GLfloat vertex[M_TNL_ATTRIB_MAX*4];	   /* current values */
   GLfloat *attrptr[M_TNL_ATTRIB_MAX];
   GLuint counter, initial_counter;
   GLboolean dangling_attr_ref;
   GLboolean have_materials;

   GLuint opcode_vertex_list;

   struct Mtnl_copied_vtx copied;

   GLfloat *current[M_TNL_ATTRIB_MAX]; /* points into ctx->ListState */
   GLubyte *currentsz[M_TNL_ATTRIB_MAX];

   void (*tabfv[M_TNL_ATTRIB_MAX][4])( const GLfloat * );
};

/** Describes an individual operation on the pipeline.
 */
struct Mtnl_pipeline_stage
{
   const char *name;
   GLuint check_state;		/* All state referenced in check() --
				 * When is the pipeline_stage struct
				 * itself invalidated?  Must be
				 * constant.
				 */

   /* Usually constant or set by the 'check' callback:
    */
   GLuint run_state;		/* All state referenced in run() --
				 * When is the cached output of the
				 * stage invalidated?
				 */

   GLboolean active;		/* True if runnable in current state */
   GLuint inputs;		/* VERT_* inputs to the stage */
   GLuint outputs;		/* VERT_* outputs of the stage */

   /* Set in _tnl_run_pipeline():
    */
   GLuint changed_inputs;	/* Generated value -- inputs to the
				 * stage that have changed since last
				 * call to 'run'.
				 */


   /* Private data for the pipeline stage:
    */
   void *privatePtr;

   /* Free private data.  May not be null.
    */
   void (*destroy)( struct Mtnl_pipeline_stage * );

   /* Called from _tnl_validate_pipeline().  Must update all fields in
    * the pipeline_stage struct for the current state.
    */
   void (*check)( MGLcontext *ctx, struct Mtnl_pipeline_stage * );

   /* Called from _tnl_run_pipeline().  The stage.changed_inputs value
    * encodes all inputs to thee struct which have changed.  If
    * non-zero, recompute all affected outputs of the stage, otherwise
    * execute any 'sideeffects' of the stage.
    *
    * Return value: GL_TRUE - keep going
    *               GL_FALSE - finished pipeline
    */
   GLboolean (*run)( MGLcontext *ctx, struct Mtnl_pipeline_stage * );
};

/** Contains the array of all pipeline stages.
 * The default values are defined at the end of t_pipeline.c */
struct Mtnl_pipeline {
   GLuint build_state_trigger;	  /**< state changes which require build */
   GLuint build_state_changes;    /**< state changes since last build */
   GLuint run_state_changes;	  /**< state changes since last run */
   GLuint run_input_changes;	  /**< VERT_* changes since last run */
   GLuint inputs;		  /**< VERT_* inputs to pipeline */
   /** This array has to end with a NULL-pointer. */
   struct Mtnl_pipeline_stage stages[MMAX_PIPELINE_STAGES+1];
   GLuint nr_stages;
};

/* Wrap all the information about vectors up in a struct.  Has
 * additional fields compared to the other vectors to help us track of
 * different vertex sizes, and whether we need to clean columns out
 * because they contain non-(0,0,0,1) values.
 *
 * The start field is used to reserve data for copied vertices at the
 * end of _mesa_transform_vb, and avoids the need for a multiplication in
 * the transformation routines.
 */
typedef struct {
   GLfloat (*data)[4];	/* may be malloc'd or point to client data */
   GLfloat *start;	/* points somewhere inside of <data> */
   GLuint count;	/* size of the vector (in elements) */
   GLuint stride;	/* stride from one element to the next (in bytes) */
   GLuint size;		/* 2-4 for vertices and 1-4 for texcoords */
   GLuint flags;	/* which columns are dirty */
   void *storage;	/* self-allocated storage */
} MGLvector4f;

/**
 * Contains the current state of a running pipeline.
 */
struct Mvertex_buffer
{
   /* Constant over life of the vertex_buffer.
    */
   GLuint      Size;

   /* Constant over the pipeline.
    */
   GLuint      Count;		              /* for everything except Elts */

   /* Pointers to current data.
    */
   GLuint      *Elts;		                
   MGLvector4f  *ObjPtr;		                /* _TNL_BIT_POS */
   MGLvector4f  *EyePtr;		                /* _TNL_BIT_POS */
   MGLvector4f  *ClipPtr;	                /* _TNL_BIT_POS */
   MGLvector4f  *NdcPtr;                         /* _TNL_BIT_POS */
   GLubyte     ClipOrMask;	                /* _TNL_BIT_POS */
   GLubyte     ClipAndMask;	                /* _TNL_BIT_POS */
   GLubyte     *ClipMask;		        /* _TNL_BIT_POS */
   MGLvector4f  *NormalPtr;	                /* _TNL_BIT_NORMAL */
   GLfloat     *NormalLengthPtr;	        /* _TNL_BIT_NORMAL */
   GLboolean   *EdgeFlag;	                /* _TNL_BIT_EDGEFLAG */
   MGLvector4f  *TexCoordPtr[MMAX_TEXTURE_COORD_UNITS]; /* VERT_TEX_0..n */
   MGLvector4f  *IndexPtr[2];	                /* _TNL_BIT_INDEX */
   MGLvector4f  *ColorPtr[2];	                /* _TNL_BIT_COLOR0 */
   MGLvector4f  *SecondaryColorPtr[2];           /* _TNL_BIT_COLOR1 */
   MGLvector4f  *PointSizePtr;	                /* _TNL_BIT_POS */
   MGLvector4f  *FogCoordPtr;	                /* _TNL_BIT_FOG */

   struct tnl_prim  *Primitive;	              
   GLuint      PrimitiveCount;	      

   /* Inputs to the vertex program stage */
   MGLvector4f *AttribPtr[M_TNL_ATTRIB_MAX];      /* GL_NV_vertex_program */

   GLuint LastClipped;
   /* Private data from _tnl_render_stage that has no business being
    * in this struct.
    */
};

struct Mtnl_vertex_arrays
{
   /* Conventional vertex attribute arrays */
   MGLvector4f  Obj;
   MGLvector4f  Normal;
   MGLvector4f  Color;
   MGLvector4f  SecondaryColor;
   MGLvector4f  FogCoord;
   MGLvector4f  TexCoord[MMAX_TEXTURE_COORD_UNITS];
   MGLvector4f  Index;

   GLubyte     *EdgeFlag;
   GLuint      *Elt;

   /* These attributes don't alias with the conventional attributes.
    * The GL_NV_vertex_program extension defines 16 extra sets of vertex
    * arrays which have precedent over the conventional arrays when enabled.
    */
   MGLvector4f  Attribs[M_TNL_ATTRIB_MAX];
};

/**
 * Context state for T&L context.
 */
typedef struct
{
   /* Driver interface.
    */
   struct Mtnl_device_driver Driver;

   /* Execute:
    */
   struct Mtnl_vtx vtx;
   
   /* Compile:
    */
   struct Mtnl_save save;

   /* Pipeline
    */
   struct Mtnl_pipeline pipeline;
   struct Mvertex_buffer vb;

   /* GLvectors for binding to vb:
    */
   struct Mtnl_vertex_arrays vtx_inputs;
   struct Mtnl_vertex_arrays save_inputs;
   struct Mtnl_vertex_arrays current;
   struct Mtnl_vertex_arrays array_inputs;

   /* Clipspace/ndc/window vertex managment:
    */
   struct Mtnl_clipspace clipspace;

   /* Probably need a better configuration mechanism:
    */
   GLboolean NeedNdcCoords;
   GLboolean LoopbackDListCassettes;
   GLboolean CalcDListNormalLengths;
   GLboolean IsolateMaterials;
   GLboolean AllowVertexFog;
   GLboolean AllowPixelFog;
   GLboolean AllowCodegen;

   GLboolean _DoVertexFog;  /* eval fog function at each vertex? */

   GLuint render_inputs;

   MGLvertexformat exec_vtxfmt;
   MGLvertexformat save_vtxfmt;

} MTNLcontext;

typedef void (GLAPIENTRY *Marray_func)( const void * );

typedef struct {
   const struct Mgl_client_array *array;
   Marray_func func;
} MAEarray;

typedef void (GLAPIENTRY *Mattrib_func)( GLuint indx, const void *data );

typedef struct {
   const struct Mgl_client_array *array;
   Mattrib_func func;
   GLuint index;
} MAEattrib;

typedef struct {
   MAEarray arrays[32];
   MAEattrib attribs[MVERT_ATTRIB_MAX + 1];
   GLuint NewState;
} MAEcontext;

#define MAE_CONTEXT(ctx) ((MAEcontext *)(ctx)->aelt_context)

#define MTNL_CONTEXT(ctx) ((MTNLcontext *)(ctx->swtnl_context))
#define MGET_VERTEX_STATE(ctx)  &(MTNL_CONTEXT(ctx)->clipspace)

#define Mforeach_s(ptr, t, list)   \
        for(ptr=(list)->next,t=(ptr)->next; list != ptr; ptr=t, t=(t)->next)
/**
 * Remove an element from list.
 *
 * \param elem element to remove.
 */
#define Mremove_from_list(elem)			\
do {						\
   (elem)->next->prev = (elem)->prev;		\
   (elem)->prev->next = (elem)->next;		\
} while (0)

/* These are used to make the ctx->Current values look like
 * arrays (with zero StrideB).
 */
struct Mac_arrays {
   struct Mgl_client_array Vertex;
   struct Mgl_client_array Normal;
   struct Mgl_client_array Color;
   struct Mgl_client_array SecondaryColor;
   struct Mgl_client_array FogCoord;
   struct Mgl_client_array Index;
   struct Mgl_client_array TexCoord[MMAX_TEXTURE_COORD_UNITS];
   struct Mgl_client_array EdgeFlag;
   struct Mgl_client_array Attrib[MVERT_ATTRIB_MAX];  /* GL_NV_vertex_program */
};

struct Mac_array_flags {
   GLboolean Vertex;
   GLboolean Normal;
   GLboolean Color;
   GLboolean SecondaryColor;
   GLboolean FogCoord;
   GLboolean Index;
   GLboolean TexCoord[MMAX_TEXTURE_COORD_UNITS];
   GLboolean EdgeFlag;
   GLboolean Attrib[MVERT_ATTRIB_MAX];  /* GL_NV_vertex_program */
};

typedef struct {
   GLuint NewState;		/* not needed? */
   GLuint NewArrayState;

   /* Facility for importing and caching array data:
    */
   struct Mac_arrays Fallback;
   struct Mac_arrays Cache;
   struct Mac_arrays Raw;
   struct Mac_array_flags IsCached;
   GLuint start;
   GLuint count;

   /* Facility for importing element lists:
    */
   GLuint *Elts;
   GLuint elt_size;

} MACcontext;

#define MAC_CONTEXT(ctx) ((MACcontext *)ctx->acache_context)

//-------------------------------------------------------------------------

/* The driver interface for the software rasterizer.
 * Unless otherwise noted, all functions are mandatory.  
 */
struct Mswrast_device_driver {

   void (*SetBuffer)(MGLcontext *ctx, MGLframebuffer *buffer, GLuint bufferBit);
   /*
    * Specifies the current color buffer for span/pixel writing/reading.
    * buffer indicates which window to write to / read from.  Normally,
    * this'll be the buffer currently bound to the context, but it doesn't
    * have to be!
    * bufferBit indicates which color buffer, exactly one of:
    *    DD_FRONT_LEFT_BIT - this buffer always exists
    *    DD_BACK_LEFT_BIT - when double buffering
    *    DD_FRONT_RIGHT_BIT - when using stereo
    *    DD_BACK_RIGHT_BIT - when using stereo and double buffering
    *    DD_AUXn_BIT - if aux buffers are implemented
    */


   /***
    *** Functions for synchronizing access to the framebuffer:
    ***/

   void (*SpanRenderStart)(MGLcontext *ctx);
   void (*SpanRenderFinish)(MGLcontext *ctx);
   /* OPTIONAL.
    *
    * Called before and after all rendering operations, including DrawPixels,
    * ReadPixels, Bitmap, span functions, and CopyTexImage, etc commands.
    * These are a suitable place for grabbing/releasing hardware locks.
    *
    * NOTE: The swrast triangle/line/point routines *DO NOT* call
    * these functions.  Locking in that case must be organized by the
    * driver by other mechanisms.
    */

   /***
    *** Functions for writing pixels to the frame buffer:
    ***/

   void (*WriteRGBASpan)( const MGLcontext *ctx,
                          GLuint n, GLint x, GLint y,
                          CONST MGLchan rgba[][4], const GLubyte mask[] );
   void (*WriteRGBSpan)( const MGLcontext *ctx,
                         GLuint n, GLint x, GLint y,
                         CONST MGLchan rgb[][3], const GLubyte mask[] );
   /* Write a horizontal run of RGBA or RGB pixels.
    * If mask is NULL, draw all pixels.
    * If mask is not null, only draw pixel [i] when mask [i] is true.
    */

   void (*WriteMonoRGBASpan)( const MGLcontext *ctx, GLuint n, GLint x, GLint y,
                              const MGLchan color[4], const GLubyte mask[] );
   /* Write a horizontal run of RGBA pixels all with the same color.
    * If mask is NULL, draw all pixels.
    * If mask is not null, only draw pixel [i] when mask [i] is true.
    */

   void (*WriteRGBAPixels)( const MGLcontext *ctx,
                            GLuint n, const GLint x[], const GLint y[],
                            CONST MGLchan rgba[][4], const GLubyte mask[] );
   /* Write array of RGBA pixels at random locations.
    */

   void (*WriteMonoRGBAPixels)( const MGLcontext *ctx,
                                GLuint n, const GLint x[], const GLint y[],
                                const MGLchan color[4], const GLubyte mask[] );
   /* Write an array of mono-RGBA pixels at random locations.
    */

   void (*WriteCI32Span)( const MGLcontext *ctx, GLuint n, GLint x, GLint y,
                          const GLuint index[], const GLubyte mask[] );
   void (*WriteCI8Span)( const MGLcontext *ctx, GLuint n, GLint x, GLint y,
                         const GLubyte index[], const GLubyte mask[] );
   /* Write a horizontal run of CI pixels.  One function is for 32bpp
    * indexes and the other for 8bpp pixels (the common case).  You mus
    * implement both for color index mode.
    * If mask is NULL, draw all pixels.
    * If mask is not null, only draw pixel [i] when mask [i] is true.
    */

   void (*WriteMonoCISpan)( const MGLcontext *ctx, GLuint n, GLint x, GLint y,
                            GLuint colorIndex, const GLubyte mask[] );
   /* Write a horizontal run of color index pixels using the color index
    * last specified by the Index() function.
    * If mask is NULL, draw all pixels.
    * If mask is not null, only draw pixel [i] when mask [i] is true.
    */

   void (*WriteCI32Pixels)( const MGLcontext *ctx,
                            GLuint n, const GLint x[], const GLint y[],
                            const GLuint index[], const GLubyte mask[] );
   /*
    * Write a random array of CI pixels.
    */

   void (*WriteMonoCIPixels)( const MGLcontext *ctx,
                              GLuint n, const GLint x[], const GLint y[],
                              GLuint colorIndex, const GLubyte mask[] );
   /* Write a random array of color index pixels using the color index
    * last specified by the Index() function.
    */


   /***
    *** Functions to read pixels from frame buffer:
    ***/

   void (*ReadCI32Span)( const MGLcontext *ctx,
                         GLuint n, GLint x, GLint y, GLuint index[] );
   /* Read a horizontal run of color index pixels.
    */

   void (*ReadRGBASpan)( const MGLcontext *ctx, GLuint n, GLint x, GLint y,
                         MGLchan rgba[][4] );
   /* Read a horizontal run of RGBA pixels.
    */

   void (*ReadCI32Pixels)( const MGLcontext *ctx,
                           GLuint n, const GLint x[], const GLint y[],
                           GLuint indx[], const GLubyte mask[] );
   /* Read a random array of CI pixels.
    */

   void (*ReadRGBAPixels)( const MGLcontext *ctx,
                           GLuint n, const GLint x[], const GLint y[],
                           MGLchan rgba[][4], const GLubyte mask[] );
   /* Read a random array of RGBA pixels.
    */



   /***
    *** For supporting hardware Z buffers:
    *** Either ALL or NONE of these functions must be implemented!
    *** NOTE that Each depth value is a 32-bit GLuint.  If the depth
    *** buffer is less than 32 bits deep then the extra upperbits are zero.
    ***/

   void (*WriteDepthSpan)( MGLcontext *ctx, GLuint n, GLint x, GLint y,
                           const MGLdepth depth[], const GLubyte mask[] );
   /* Write a horizontal span of values into the depth buffer.  Only write
    * depth[i] value if mask[i] is nonzero.
    */

   void (*ReadDepthSpan)( MGLcontext *ctx, GLuint n, GLint x, GLint y,
                          MGLdepth depth[] );
   /* Read a horizontal span of values from the depth buffer.
    */


   void (*WriteDepthPixels)( MGLcontext *ctx, GLuint n,
                             const GLint x[], const GLint y[],
                             const MGLdepth depth[], const GLubyte mask[] );
   /* Write an array of randomly positioned depth values into the
    * depth buffer.  Only write depth[i] value if mask[i] is nonzero.
    */

   void (*ReadDepthPixels)( MGLcontext *ctx, GLuint n,
                            const GLint x[], const GLint y[],
                            MGLdepth depth[] );
   /* Read an array of randomly positioned depth values from the depth buffer.
    */



   /***
    *** For supporting hardware stencil buffers:
    *** Either ALL or NONE of these functions must be implemented!
    ***/

   void (*WriteStencilSpan)( MGLcontext *ctx, GLuint n, GLint x, GLint y,
                             const MGLstencil stencil[], const GLubyte mask[] );
   /* Write a horizontal span of stencil values into the stencil buffer.
    * If mask is NULL, write all stencil values.
    * Else, only write stencil[i] if mask[i] is non-zero.
    */

   void (*ReadStencilSpan)( MGLcontext *ctx, GLuint n, GLint x, GLint y,
                            MGLstencil stencil[] );
   /* Read a horizontal span of stencil values from the stencil buffer.
    */

   void (*WriteStencilPixels)( MGLcontext *ctx, GLuint n,
                               const GLint x[], const GLint y[],
                               const MGLstencil stencil[],
                               const GLubyte mask[] );
   /* Write an array of stencil values into the stencil buffer.
    * If mask is NULL, write all stencil values.
    * Else, only write stencil[i] if mask[i] is non-zero.
    */

   void (*ReadStencilPixels)( MGLcontext *ctx, GLuint n,
                              const GLint x[], const GLint y[],
                              MGLstencil stencil[] );
   /* Read an array of stencil values from the stencil buffer.
    */
};

/**
 * \struct SWvertex
 * \brief Data-structure to handle vertices in the software rasterizer.
 * 
 * The software rasterizer now uses this format for vertices.  Thus a
 * 'RasterSetup' stage or other translation is required between the
 * tnl module and the swrast rasterization functions.  This serves to
 * isolate the swrast module from the internals of the tnl module, and
 * improve its usefulness as a fallback mechanism for hardware
 * drivers.
 *
 * Full software drivers:
 *   - Register the rastersetup and triangle functions from
 *     utils/software_helper.
 *   - On statechange, update the rasterization pointers in that module.
 *
 * Rasterization hardware drivers:
 *   - Keep native rastersetup.
 *   - Implement native twoside,offset and unfilled triangle setup.
 *   - Implement a translator from native vertices to swrast vertices.
 *   - On partial fallback (mix of accelerated and unaccelerated
 *   prims), call a pass-through function which translates native
 *   vertices to SWvertices and calls the appropriate swrast function.
 *   - On total fallback (vertex format insufficient for state or all
 *     primitives unaccelerated), hook in swrast_setup instead.
 */
typedef struct {
   /** win[0], win[1] are the screen-coords of SWvertex. win[2] is the
    * z-coord. what is win[3]? */
   GLfloat win[4];
   GLfloat texcoord[MMAX_TEXTURE_COORD_UNITS][4];
   MGLchan color[4];
   MGLchan specular[4];
   GLfloat fog;
   GLfloat index;
   GLfloat pointSize;
} MSWvertex;

typedef void (*Mswrast_point_func)( MGLcontext *ctx, const MSWvertex *);
typedef void (*Mswrast_line_func)( MGLcontext *ctx, const MSWvertex *, const MSWvertex *);
typedef void (*Mswrast_tri_func)( MGLcontext *ctx, const MSWvertex *, const MSWvertex *, const MSWvertex *);
typedef void (_ASMAPIP Mblend_func)( MGLcontext *ctx, GLuint n, const GLubyte mask[], MGLchan src[][4], CONST MGLchan dst[][4] );
typedef void (*Mtexture_sample_func)(MGLcontext *ctx, GLuint texUnit, const struct Mgl_texture_object *tObj, GLuint n, const GLfloat texcoords[][4], const GLfloat lambda[], MGLchan rgba[][4]);

/**
 * \struct sw_span
 * \brief Contains data for either a horizontal line or a set of
 * pixels that are passed through a pipeline of functions before being
 * drawn.
 *
 * The sw_span structure describes the colors, Z, fogcoord, texcoords,
 * etc for either a horizontal run or an array of independent pixels.
 * We can either specify a base/step to indicate interpolated values, or
 * fill in arrays of values.  The interpMask and arrayMask bitfields
 * indicate which are active.
 *
 * With this structure it's easy to hand-off span rasterization to
 * subroutines instead of doing it all inline in the triangle functions
 * like we used to do.
 * It also cleans up the local variable namespace a great deal.
 *
 * It would be interesting to experiment with multiprocessor rasterization
 * with this structure.  The triangle rasterizer could simply emit a
 * stream of these structures which would be consumed by one or more
 * span-processing threads which could run in parallel.
 */
struct Msw_span {
   GLint x, y;

   /** Only need to process pixels between start <= i < end */
   /** At this time, start is always zero. */
   GLuint start, end;

   /** This flag indicates that mask[] array is effectively filled with ones */
   GLboolean writeAll;

   /** either GL_POLYGON, GL_LINE, GL_POLYGON, GL_BITMAP */
   GLenum primitive;

   /** 0 = front-facing span, 1 = back-facing span (for two-sided stencil) */
   GLuint facing;

   /**
    * This bitmask (of  \link SpanFlags SPAN_* flags\endlink) indicates
    * which of the x/xStep variables are relevant.
    */
   GLuint interpMask;

   /* For horizontal spans, step is the partial derivative wrt X.
    * For lines, step is the delta from one fragment to the next.
    */
#if CHAN_TYPE == GL_FLOAT
   GLfloat red, redStep;
   GLfloat green, greenStep;
   GLfloat blue, blueStep;
   GLfloat alpha, alphaStep;
   GLfloat specRed, specRedStep;
   GLfloat specGreen, specGreenStep;
   GLfloat specBlue, specBlueStep;
#else /* CHAN_TYPE == GL_UNSIGNED_BYTE or GL_UNSIGNED_SHORT */
   MGLfixed red, redStep;
   MGLfixed green, greenStep;
   MGLfixed blue, blueStep;
   MGLfixed alpha, alphaStep;
   MGLfixed specRed, specRedStep;
   MGLfixed specGreen, specGreenStep;
   MGLfixed specBlue, specBlueStep;
#endif
   MGLfixed index, indexStep;
   MGLfixed z, zStep;
   GLfloat fog, fogStep;
   GLfloat tex[MMAX_TEXTURE_COORD_UNITS][4];  /* s, t, r, q */
   GLfloat texStepX[MMAX_TEXTURE_COORD_UNITS][4];
   GLfloat texStepY[MMAX_TEXTURE_COORD_UNITS][4];
   MGLfixed intTex[2], intTexStep[2];  /* s, t only */

   /* partial derivatives wrt X and Y. */
   GLfloat dzdx, dzdy;
   GLfloat w, dwdx, dwdy;
   GLfloat drdx, drdy;
   GLfloat dgdx, dgdy;
   GLfloat dbdx, dbdy;
   GLfloat dadx, dady;
   GLfloat dsrdx, dsrdy;
   GLfloat dsgdx, dsgdy;
   GLfloat dsbdx, dsbdy;
   GLfloat dfogdx, dfogdy;

   /**
    * This bitmask (of \link SpanFlags SPAN_* flags\endlink) indicates
    * which of the fragment arrays in the span_arrays struct are relevant.
    */
   GLuint arrayMask;

   /**
    * We store the arrays of fragment values in a separate struct so
    * that we can allocate sw_span structs on the stack without using
    * a lot of memory.  The span_arrays struct is about 400KB while the
    * sw_span struct is only about 512 bytes.
    */
   struct span_arrays *array;
};

/**
 * \struct SWcontext
 * \brief SWContext?
 */
typedef struct
{
   /** Driver interface:
    */
   struct Mswrast_device_driver Driver;

   /** Configuration mechanisms to make software rasterizer match
    * characteristics of the hardware rasterizer (if present):
    */
   GLboolean AllowVertexFog;
   GLboolean AllowPixelFog;

   /** Derived values, invalidated on statechanges, updated from
    * _swrast_validate_derived():
    */
   GLuint _RasterMask;
   GLfloat _MinMagThresh[MMAX_TEXTURE_IMAGE_UNITS];
   GLfloat _BackfaceSign;
   GLboolean _PreferPixelFog;    /* Compute fog blend factor per fragment? */
   GLboolean _AnyTextureCombine;
   MGLchan _FogColor[3];
   GLboolean _FogEnabled;

   /* Accum buffer temporaries.
    */
   GLboolean _IntegerAccumMode;	/**< Storing unscaled integers? */
   GLfloat _IntegerAccumScaler;	/**< Implicit scale factor */

   MGLchan *CurAuxBuffer;

   /* Working values:
    */
   GLuint StippleCounter;    /**< Line stipple counter */
   GLuint NewState;
   GLuint StateChanges;
   GLenum Primitive;    /* current primitive being drawn (ala glBegin) */
   GLbitfield CurrentBufferBit; /* exactly one the of DD_*_BIT buffer bits */

   /** Mechanism to allow driver (like X11) to register further
    * software rasterization routines.
    */
   /*@{*/
   void (*choose_point)( MGLcontext * );
   void (*choose_line)( MGLcontext * );
   void (*choose_triangle)( MGLcontext * );

   GLuint invalidate_point;
   GLuint invalidate_line;
   GLuint invalidate_triangle;
   /*@}*/

   /** Function pointers for dispatch behind public entrypoints. */
   /*@{*/
   void (*InvalidateState)( MGLcontext *ctx, GLuint new_state );

   Mswrast_point_func Point;
   Mswrast_line_func Line;
   Mswrast_tri_func Triangle;
   /*@}*/

   /**
    * Placeholders for when separate specular (or secondary color) is
    * enabled but texturing is not.
    */
   /*@{*/
   Mswrast_point_func SpecPoint;
   Mswrast_line_func SpecLine;
   Mswrast_tri_func SpecTriangle;
   /*@}*/

   /**
    * Typically, we'll allocate a sw_span structure as a local variable
    * and set its 'array' pointer to point to this object.  The reason is
    * this object is big and causes problems when allocated on the stack
    * on some systems.
    */
   struct Mspan_arrays *SpanArrays;

   /**
    * Used to buffer N GL_POINTS, instead of rendering one by one.
    */
   struct Msw_span PointSpan;

   /** Internal hooks, kept uptodate by the same mechanism as above.
    */
   Mblend_func BlendFunc;
   Mtexture_sample_func TextureSample[MMAX_TEXTURE_IMAGE_UNITS];

   /** Buffer for saving the sampled texture colors.
    * Needed for GL_ARB_texture_env_crossbar implementation.
    */
   MGLchan *TexelBuffer;

} MSWcontext;

//extern void _swrast_validate_derived( GLcontext *ctx );
#define MSWRAST_CONTEXT(ctx) ((MSWcontext *)ctx->swrast_context)

/**
 * Frame buffer.
 *
 * A "frame buffer" is a color buffer and its optional ancillary buffers:
 * depth, accum, stencil, and software-simulated alpha buffers.
 * In C++ terms, think of this as a base class from which device drivers
 * will make derived classes.
 */
struct Mgl_frame_buffer
{
   MGLvisual Visual;		/**< The corresponding visual */

   GLuint Width, Height;	/**< size of frame buffer in pixels */

   GLboolean UseSoftwareDepthBuffer;
   GLboolean UseSoftwareAccumBuffer;
   GLboolean UseSoftwareStencilBuffer;
   GLboolean UseSoftwareAlphaBuffers;
   GLboolean UseSoftwareAuxBuffers;

   /** \name Software depth (aka Z) buffer */
   /*@{*/
   GLvoid *DepthBuffer;		/**< array [Width*Height] of GLushort or GLuint*/
   /*@}*/

   /** \name Software stencil buffer */
   /*@{*/
   MGLstencil *Stencil;		/**< array [Width*Height] of GLstencil values */
   /*@}*/

   /** \name Software accumulation buffer */
   /*@{*/
   MGLaccum *Accum;		/**< array [4*Width*Height] of GLaccum values */
   /*@}*/

   /** \name Software alpha planes */
   /*@{*/
   MGLchan *FrontLeftAlpha;	/**< array [Width*Height] of GLchan */
   MGLchan *BackLeftAlpha;	/**< array [Width*Height] of GLchan */
   MGLchan *FrontRightAlpha;	/**< array [Width*Height] of GLchan */
   MGLchan *BackRightAlpha;	/**< array [Width*Height] of GLchan */
   /*@}*/

   MGLchan *AuxBuffers[MMAX_AUX_BUFFERS];

   /** 
    * \name Drawing bounds
    *
    * Intersection of window size and scissor box 
    */
   /*@{*/
   GLint _Xmin;  /**< inclusive */
   GLint _Ymin;  /**< inclusive */
   GLint _Xmax;  /**< exclusive */
   GLint _Ymax;  /**< exclusive */
   /*@}*/
};

//-------------------------------------------------------------------------

extern void MOSMesaDestroyContext(MOSMesaContext ctx);
extern MOSMesaContext MOSMesaCreateContextExt(GLenum format, GLint depthBits, GLint stencilBits, GLint accumBits, MOSMesaContext sharelist);
extern GLboolean MOSMesaMakeCurrent(MOSMesaContext ctx, void *buffer, GLenum type, GLsizei width, GLsizei height);

//-------------------------------------------------------------------------

#endif