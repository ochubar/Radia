
#ifndef MOSMESA_H
#include "mosmesa.h"
#endif

#include <stdlib.h>

//-------------------------------------------------------------------------
/** Wrapper around either free() or xf86free() */
void M_mesa_free(void *ptr)
{
#if defined(XFree86LOADER) && defined(IN_MODULE)
   xf86free(ptr);
#else
   free(ptr);
#endif
}

/**
 * Free memory allocated with _mesa_align_malloc() or _mesa_align_calloc().
 *
 * \param ptr pointer to the memory to be freed.
 * 
 * The actual address to free is stored in the word immediately before the
 * address the client sees.
 */
void M_mesa_align_free(void *ptr)
{
#if 0
   M_mesa_free( (void *)(*(unsigned long *)((unsigned long)ptr - sizeof(void *))) );
#else
   void **cubbyHole = (void **) ((char *) ptr - sizeof(void *));
   void *realAddr = *cubbyHole;
   M_mesa_free(realAddr);
#endif
}

/** Free memory */
#define MFREE(PTR)          M_mesa_free(PTR)
#define MALIGN_FREE(PTR)            M_mesa_align_free(PTR)

#ifdef MESA_EXTERNAL_BUFFERALLOC
extern void *M_ext_mesa_alloc_pixelbuffer( unsigned int size );
extern void M_ext_mesa_free_pixelbuffer( void *pb );

#define MESA_PBUFFER_ALLOC(BYTES)  (void *) M_ext_mesa_alloc_pixelbuffer(BYTES)
#define MESA_PBUFFER_FREE(PTR)     M_ext_mesa_free_pixelbuffer(PTR)
#else
/* Default buffer allocation uses the aligned allocation routines: */
#define MESA_PBUFFER_ALLOC(BYTES)  (void *) M_mesa_align_malloc(BYTES, 512)
#define MESA_PBUFFER_FREE(PTR)     M_mesa_align_free(PTR)
#endif

//-------------------------------------------------------------------------

/* Used when thread safety disabled */
void *M_glapi_Context = NULL;

//-------------------------------------------------------------------------

void M_tnl_free_vertices( MGLcontext *ctx )
{
   struct Mtnl_clipspace *vtx = MGET_VERTEX_STATE(ctx);

   if (vtx->vertex_buf) {
      MALIGN_FREE(vtx->vertex_buf);
      vtx->vertex_buf = 0;
   }
}

//-------------------------------------------------------------------------

void M_swsetup_DestroyContext(MGLcontext *ctx)
{
   MSScontext *swsetup = MSWSETUP_CONTEXT(ctx);

   if (swsetup) {
      MFREE(swsetup);
      ctx->swsetup_context = 0;
   }

   M_tnl_free_vertices( ctx );
}

//-------------------------------------------------------------------------
/**
 * Destroy the context's vertex array stuff.
 * Called during T 'n L context destruction.
 */
void M_tnl_array_destroy( MGLcontext *ctx )
{
}

//-------------------------------------------------------------------------
/**
 * Deallocate the immediate-mode buffer for the given context, if
 * its reference count goes to zero.
 */
void M_tnl_save_destroy( MGLcontext *ctx )
{
}

//-------------------------------------------------------------------------

static void Mfree_funcs( struct M_tnl_dynfn *l )
{
   struct M_tnl_dynfn *f, *tmp;
   Mforeach_s (f, tmp, l) {
      Mremove_from_list( f );
      MALIGN_FREE( f->code );
      MFREE( f );
   }
}

//-------------------------------------------------------------------------

void M_tnl_vtx_destroy( MGLcontext *ctx )
{
   MTNLcontext *tnl = MTNL_CONTEXT(ctx); 
   GLuint i;

   for (i = 0; i < 4; i++) {
      Mfree_funcs( &tnl->vtx.cache.Vertex[i] );
      Mfree_funcs( &tnl->vtx.cache.Attribute[i] ); 
   }
}

//-------------------------------------------------------------------------

void M_tnl_destroy_pipeline( MGLcontext *ctx )
{
	MTNLcontext *tnl = MTNL_CONTEXT(ctx);
	GLuint i;

	for (i = 0 ; i < tnl->pipeline.nr_stages ; i++)
		tnl->pipeline.stages[i].destroy( &tnl->pipeline.stages[i] );

	tnl->pipeline.nr_stages = 0;
}

//-------------------------------------------------------------------------

void M_ae_destroy_context( MGLcontext *ctx )
{
	if ( MAE_CONTEXT( ctx ) ) {
		MFREE( ctx->aelt_context );
		ctx->aelt_context = 0;
	}
}

//-------------------------------------------------------------------------

void M_tnl_DestroyContext( MGLcontext *ctx )
{
    MTNLcontext *tnl = MTNL_CONTEXT(ctx);
    M_tnl_array_destroy( ctx );
    M_tnl_vtx_destroy( ctx );
	M_tnl_save_destroy( ctx );
	M_tnl_destroy_pipeline( ctx );
    M_ae_destroy_context( ctx );
    MFREE(tnl);
    ctx->swtnl_context = 0;
}

//-------------------------------------------------------------------------

void M_ac_DestroyContext( MGLcontext *ctx )
{
	struct Mgl_buffer_object *nullObj = ctx->Array.NullBufferObj;
	MACcontext *ac = MAC_CONTEXT(ctx);
	GLint i;

	//only free vertex data if it's really a pointer to vertex data and not an offset into a buffer object.
	if(ac->Cache.Vertex.Ptr && ac->Cache.Vertex.BufferObj == nullObj) MFREE( (void *) ac->Cache.Vertex.Ptr );
	if(ac->Cache.Normal.Ptr && ac->Cache.Normal.BufferObj == nullObj) MFREE( (void *) ac->Cache.Normal.Ptr );
	if (ac->Cache.Color.Ptr && ac->Cache.Color.BufferObj == nullObj)
		MFREE( (void *) ac->Cache.Color.Ptr );
	if (ac->Cache.SecondaryColor.Ptr && ac->Cache.SecondaryColor.BufferObj == nullObj)
		MFREE( (void *) ac->Cache.SecondaryColor.Ptr );
	if (ac->Cache.EdgeFlag.Ptr && ac->Cache.EdgeFlag.BufferObj == nullObj)
		MFREE( (void *) ac->Cache.EdgeFlag.Ptr );
	if (ac->Cache.Index.Ptr && ac->Cache.Index.BufferObj == nullObj)
		MFREE( (void *) ac->Cache.Index.Ptr );
	if (ac->Cache.FogCoord.Ptr && ac->Cache.FogCoord.BufferObj == nullObj)
		MFREE( (void *) ac->Cache.FogCoord.Ptr );

	for (i = 0; i < MMAX_TEXTURE_COORD_UNITS; i++) {
		if (ac->Cache.TexCoord[i].Ptr && ac->Cache.TexCoord[i].BufferObj == nullObj)
			MFREE( (void *) ac->Cache.TexCoord[i].Ptr );
	}

	for (i = 0; i < MVERT_ATTRIB_MAX; i++) {
		if (ac->Cache.Attrib[i].Ptr && ac->Cache.Attrib[i].BufferObj == nullObj)
			MFREE( (void *) ac->Cache.Attrib[i].Ptr );
	}

    if(ac->Elts) MFREE( ac->Elts );

   /* Free the context structure itself */
    MFREE(ac);
    ctx->acache_context = NULL;
}

//-------------------------------------------------------------------------

void M_swrast_DestroyContext( MGLcontext *ctx )
{
	MSWcontext *swrast = MSWRAST_CONTEXT(ctx);

	//if (SWRAST_DEBUG) {
	//	_mesa_debug(ctx, "_swrast_DestroyContext\n");
	//}

	MFREE( swrast->SpanArrays );
	MFREE( swrast->TexelBuffer );
	MFREE( swrast );
	ctx->swrast_context = 0;
}

//-------------------------------------------------------------------------
/**
 * Destroy a visual and free its memory.
 *
 * \param vis visual.
 * 
 * Frees the visual structure.
 */
void M_mesa_destroy_visual( MGLvisual *vis )
{
    MFREE(vis);
}

//-------------------------------------------------------------------------
/**
 * Free the data hanging off of \p buffer, but not \p buffer itself.
 *
 * \param buffer framebuffer.
 *
 * Frees all the buffers associated with the structure.
 */
void M_mesa_free_framebuffer_data( MGLframebuffer *buffer )
{
   if (!buffer)
      return;

   if (buffer->UseSoftwareDepthBuffer && buffer->DepthBuffer) {
      MESA_PBUFFER_FREE( buffer->DepthBuffer );
      buffer->DepthBuffer = NULL;
   }
   if (buffer->UseSoftwareAccumBuffer && buffer->Accum) {
      MESA_PBUFFER_FREE( buffer->Accum );
      buffer->Accum = NULL;
   }
   if (buffer->UseSoftwareStencilBuffer && buffer->Stencil) {
      MESA_PBUFFER_FREE( buffer->Stencil );
      buffer->Stencil = NULL;
   }
   if (buffer->UseSoftwareAlphaBuffers){
      if (buffer->FrontLeftAlpha) {
         MESA_PBUFFER_FREE( buffer->FrontLeftAlpha );
         buffer->FrontLeftAlpha = NULL;
      }
      if (buffer->BackLeftAlpha) {
         MESA_PBUFFER_FREE( buffer->BackLeftAlpha );
         buffer->BackLeftAlpha = NULL;
      }
      if (buffer->FrontRightAlpha) {
         MESA_PBUFFER_FREE( buffer->FrontRightAlpha );
         buffer->FrontRightAlpha = NULL;
      }
      if (buffer->BackRightAlpha) {
         MESA_PBUFFER_FREE( buffer->BackRightAlpha );
         buffer->BackRightAlpha = NULL;
      }
   }
}

//-------------------------------------------------------------------------
/**
 * Free a framebuffer struct and its buffers.
 *
 * Calls _mesa_free_framebuffer_data() and frees the structure.
 */
void M_mesa_destroy_framebuffer( MGLframebuffer *buffer )
{
   if (buffer) {
      M_mesa_free_framebuffer_data(buffer);
      MFREE(buffer);
   }
}

//-------------------------------------------------------------------------
/*
 * Get the current context pointer for this thread.
 * The context pointer is an opaque type which should be cast from
 * void to the real context pointer type.
 */
void *M_glapi_get_context(void)
{
#if defined(THREADS)
   if (ThreadSafe) {
      return _glthread_GetTSD(&ContextTSD);
   }
   else {
      return M_glapi_Context;
   }
#else
   return M_glapi_Context;
#endif
}

//-------------------------------------------------------------------------
/**
 * Get current context for the calling thread.
 * 
 * \return pointer to the current GL context.
 * 
 * Calls _glapi_get_context(). This isn't the fastest way to get the current
 * context.  If you need speed, see the #GET_CURRENT_CONTEXT macro in context.h.
 */
MGLcontext *M_mesa_get_current_context( void )
{
   return (MGLcontext *) M_glapi_get_context();
}

//-------------------------------------------------------------------------
/**
 * Check if the given context can render into the given framebuffer
 * by checking visual attributes.
 * \return GL_TRUE if compatible, GL_FALSE otherwise.
 */
static GLboolean Mcheck_compatible(const MGLcontext *ctx, const MGLframebuffer *buffer)
{
   const MGLvisual *ctxvis = &ctx->Visual;
   const MGLvisual *bufvis = &buffer->Visual;

   if (ctxvis == bufvis)
      return GL_TRUE;

   if (ctxvis->rgbMode != bufvis->rgbMode)
      return GL_FALSE;
   if (ctxvis->doubleBufferMode && !bufvis->doubleBufferMode)
      return GL_FALSE;
   if (ctxvis->stereoMode && !bufvis->stereoMode)
      return GL_FALSE;
   if (ctxvis->haveAccumBuffer && !bufvis->haveAccumBuffer)
      return GL_FALSE;
   if (ctxvis->haveDepthBuffer && !bufvis->haveDepthBuffer)
      return GL_FALSE;
   if (ctxvis->haveStencilBuffer && !bufvis->haveStencilBuffer)
      return GL_FALSE;
   if (ctxvis->redMask && ctxvis->redMask != bufvis->redMask)
      return GL_FALSE;
   if (ctxvis->greenMask && ctxvis->greenMask != bufvis->greenMask)
      return GL_FALSE;
   if (ctxvis->blueMask && ctxvis->blueMask != bufvis->blueMask)
      return GL_FALSE;
   if (ctxvis->depthBits && ctxvis->depthBits != bufvis->depthBits)
      return GL_FALSE;
   if (ctxvis->stencilBits && ctxvis->stencilBits != bufvis->stencilBits)
      return GL_FALSE;

   return GL_TRUE;
}

//-------------------------------------------------------------------------
/*
 * Set the global or per-thread dispatch table pointer.
 */
//void M_glapi_set_dispatch(struct M_glapi_table *dispatch)
//{
//   if (!dispatch) {
//      /* use the no-op functions */
//      dispatch = (struct M_glapi_table *) __glapi_noop_table;
//   }
//#ifdef DEBUG
//   else {
//      _glapi_check_table(dispatch);
//   }
//#endif
//
//#if defined(THREADS)
//   if (DispatchOverride) {
//      _glthread_SetTSD(&RealDispatchTSD, (void *) dispatch);
//      if (ThreadSafe)
//         _glapi_RealDispatch = (struct _glapi_table*) __glapi_threadsafe_table;
//      else
//         _glapi_RealDispatch = dispatch;
//   }
//   else {
//      /* normal operation */
//      _glthread_SetTSD(&_gl_DispatchTSD, (void *) dispatch);
//      if (ThreadSafe) {
//	 _glapi_Dispatch = (struct _glapi_table *) __glapi_threadsafe_table;
//	 _glapi_DispatchTSD = NULL;
//      }
//      else {
//	 _glapi_Dispatch = dispatch;
//	 _glapi_DispatchTSD = dispatch;
//      }
//   }
//#else /*THREADS*/
//   if (DispatchOverride) {
//      _glapi_RealDispatch = dispatch;
//   }
//   else {
//      _glapi_Dispatch = dispatch;
//   }
//#endif /*THREADS*/
//}

//-------------------------------------------------------------------------
/*
 * We should call this periodically from a function such as glXMakeCurrent
 * in order to test if multiple threads are being used.
 */
void M_glapi_check_multithread(void)
{
//#if defined(THREADS)
   if (!ThreadSafe) {
      static unsigned long knownID;
      static GLboolean firstCall = GL_TRUE;
      if (firstCall) {
         knownID = _glthread_GetID();
         firstCall = GL_FALSE;
      }
      else if (knownID != _glthread_GetID()) {
         ThreadSafe = GL_TRUE;
         _glapi_set_dispatch(NULL);
      }
   }
   else if (!_glapi_get_dispatch()) {
      /* make sure that this thread's dispatch pointer isn't null */
      _glapi_set_dispatch(NULL);
   }
//#endif
}

//-------------------------------------------------------------------------
/**
 * Bind the given context to the given draw-buffer and read-buffer and
 * make it the current context for this thread.
 *
 * \param newCtx new GL context. If NULL then there will be no current GL
 * context.
 * \param drawBuffer draw framebuffer.
 * \param readBuffer read framebuffer.
 * 
 * Check that the context's and framebuffer's visuals are compatible, returning
 * immediately otherwise. Sets the glapi current context via
 * _glapi_set_context(). If \p newCtx is not NULL, associates \p drawBuffer and
 * \p readBuffer with it and calls dd_function_table::ResizeBuffers if the buffers size has changed. 
 * Calls dd_function_table::MakeCurrent callback if defined.
 *
 * When a context is bound by the first time and the \c MESA_INFO environment
 * variable is set it calls print_info() as an aid for remote user
 * troubleshooting.
 */
void M_mesa_make_current2( MGLcontext *newCtx, MGLframebuffer *drawBuffer, MGLframebuffer *readBuffer )
{
   //if (MESA_VERBOSE)
   //   _mesa_debug(newCtx, "_mesa_make_current2()\n");

   // Check that the context's and framebuffer's visuals are compatible.
   if (newCtx && drawBuffer && newCtx->DrawBuffer != drawBuffer) {
      if (!Mcheck_compatible(newCtx, drawBuffer))
         return;
   }
   if (newCtx && readBuffer && newCtx->ReadBuffer != readBuffer) {
      if (!Mcheck_compatible(newCtx, readBuffer))
         return;
   }

   // We call this function periodically (just here for now) in
   // order to detect when multithreading has begun.
   M_glapi_check_multithread();

/*

   _glapi_set_context((void *) newCtx);
   ASSERT(_mesa_get_current_context() == newCtx);


   if (!newCtx) {
      _glapi_set_dispatch(NULL);  // none current
   }
   else {
      _glapi_set_dispatch(newCtx->CurrentDispatch);

      if (drawBuffer && readBuffer) {
	 // TODO: check if newCtx and buffer's visual match???
	 newCtx->DrawBuffer = drawBuffer;
	 newCtx->ReadBuffer = readBuffer;
	 newCtx->NewState |= _NEW_BUFFERS;

#if _HAVE_FULL_GL
         if (drawBuffer->Width == 0 && drawBuffer->Height == 0) {
            // get initial window size
            GLuint bufWidth, bufHeight;

            // ask device driver for size of output buffer
            (*newCtx->Driver.GetBufferSize)( drawBuffer, &bufWidth, &bufHeight );

            if (drawBuffer->Width != bufWidth || drawBuffer->Height != bufHeight) {

	       drawBuffer->Width = bufWidth;
	       drawBuffer->Height = bufHeight;

 	       newCtx->Driver.ResizeBuffers( drawBuffer );
	    }
         }

         if (readBuffer != drawBuffer &&
             readBuffer->Width == 0 && readBuffer->Height == 0) {
            // get initial window size 
            GLuint bufWidth, bufHeight;

            // ask device driver for size of output buffer
            (*newCtx->Driver.GetBufferSize)( readBuffer, &bufWidth, &bufHeight );

            if (readBuffer->Width != bufWidth ||
		readBuffer->Height != bufHeight) {

	       readBuffer->Width = bufWidth;
	       readBuffer->Height = bufHeight;

	       newCtx->Driver.ResizeBuffers( readBuffer );
	    }
         }
#endif
      }

      // Alert the driver - usually passed on to the sw t&l module,
      // but also used to detect threaded cases in the radeon codegen
      // hw t&l module.
      if (newCtx->Driver.MakeCurrent)
	 newCtx->Driver.MakeCurrent( newCtx, drawBuffer, readBuffer );

      // We can use this to help debug user's problems.  Tell them to set
      // the MESA_INFO env variable before running their app.  Then the
      // first time each context is made current we'll print some useful
      // information.
      if (newCtx->FirstTimeCurrent) {
	 if (_mesa_getenv("MESA_INFO")) {
	    _mesa_print_info();
	 }
	 newCtx->FirstTimeCurrent = GL_FALSE;
      }
   }
*/
}

/**
 * Set the current context, binding the given frame buffer to the context.
 *
 * \param newCtx new GL context.
 * \param buffer framebuffer.
 * 
 * Calls _mesa_make_current2() with \p buffer as read and write framebuffer.
 */
void M_mesa_make_current( MGLcontext *newCtx, MGLframebuffer *buffer )
{
	M_mesa_make_current2( newCtx, buffer, buffer );
}

//-------------------------------------------------------------------------
/**
 * Free the data associated with the given context.
 * 
 * But doesn't free the GLcontext struct itself.
 *
 * \sa _mesa_initialize_context() and init_attrib_groups().
 */
void M_mesa_free_context_data( MGLcontext *ctx )
{
   // if we're destroying the current context, unbind it first
   if (ctx == M_mesa_get_current_context()) {
      M_mesa_make_current(NULL, NULL);
   }
/*
   _mesa_free_lighting_data( ctx );
   _mesa_free_eval_data( ctx );
   _mesa_free_texture_data( ctx );
   _mesa_free_matrix_data( ctx );
   _mesa_free_viewport_data( ctx );
   _mesa_free_colortables_data( ctx );
   _mesa_free_program_data(ctx);
   _mesa_free_occlude_data(ctx);

#if FEATURE_ARB_vertex_buffer_object
   _mesa_delete_buffer_object(ctx, ctx->Array.NullBufferObj);
#endif

   // Shared context state (display lists, textures, etc)
   _glthread_LOCK_MUTEX(ctx->Shared->Mutex);
   ctx->Shared->RefCount--;
   assert(ctx->Shared->RefCount >= 0);
   _glthread_UNLOCK_MUTEX(ctx->Shared->Mutex);
   if (ctx->Shared->RefCount == 0) {
      // free shared state
      free_shared_state( ctx, ctx->Shared );
   }
*/
   if (ctx->Extensions.String) MFREE((void *) ctx->Extensions.String);

   MFREE(ctx->Exec);
   MFREE(ctx->Save);
}

//-------------------------------------------------------------------------
/*
 * Destroy an Off-Screen Mesa rendering context.
 *
 * Input:  ctx - the context to destroy
 */
void MOSMesaDestroyContext(MOSMesaContext ctx)
{
    if(ctx) 
	{
        M_swsetup_DestroyContext( &ctx->mesa );
        M_tnl_DestroyContext( &ctx->mesa );
		M_ac_DestroyContext( &ctx->mesa );
		M_swrast_DestroyContext( &ctx->mesa );
		M_mesa_destroy_visual( ctx->gl_visual );
        M_mesa_destroy_framebuffer( ctx->gl_buffer );
		M_mesa_free_context_data( &ctx->mesa );

		MFREE( ctx );
	}
}

//-------------------------------------------------------------------------
/*
 * New in Mesa 3.5
 *
 * Create context and specify size of ancillary buffers.
 */
MOSMesaContext MOSMesaCreateContextExt(GLenum format, GLint depthBits, GLint stencilBits, GLint accumBits, MOSMesaContext sharelist)
{
	MOSMesaContext osmesa;

	struct Mdd_function_table functions;
	GLint rshift, gshift, bshift, ashift;
	GLint rind, gind, bind, aind;
	GLint indexBits = 0, redBits = 0, greenBits = 0, blueBits = 0, alphaBits =0;
	GLboolean rgbmode;
	const GLuint i4 = 1;
	const GLubyte *i1 = (GLubyte *) &i4;
	const GLint little_endian = *i1;

    rind = gind = bind = aind = 0;
/*
   if (format==OSMESA_COLOR_INDEX) {
      indexBits = 8;
      rshift = gshift = bshift = ashift = 0;
      rgbmode = GL_FALSE;
   }
   else if (format==OSMESA_RGBA) {
      indexBits = 0;
      redBits = CHAN_BITS;
      greenBits = CHAN_BITS;
      blueBits = CHAN_BITS;
      alphaBits = CHAN_BITS;
      rind = 0;
      gind = 1;
      bind = 2;
      aind = 3;
      if (little_endian) {
         rshift = 0;
         gshift = 8;
         bshift = 16;
         ashift = 24;
      }
      else {
         rshift = 24;
         gshift = 16;
         bshift = 8;
         ashift = 0;
      }
      rgbmode = GL_TRUE;
   }
   else if (format==OSMESA_BGRA) {
      indexBits = 0;
      redBits = CHAN_BITS;
      greenBits = CHAN_BITS;
      blueBits = CHAN_BITS;
      alphaBits = CHAN_BITS;
      bind = 0;
      gind = 1;
      rind = 2;
      aind = 3;
      if (little_endian) {
         bshift = 0;
         gshift = 8;
         rshift = 16;
         ashift = 24;
      }
      else {
         bshift = 24;
         gshift = 16;
         rshift = 8;
         ashift = 0;
      }
      rgbmode = GL_TRUE;
   }
   else if (format==OSMESA_ARGB) {
      indexBits = 0;
      redBits = CHAN_BITS;
      greenBits = CHAN_BITS;
      blueBits = CHAN_BITS;
      alphaBits = CHAN_BITS;
      aind = 0;
      rind = 1;
      gind = 2;
      bind = 3;
      if (little_endian) {
         ashift = 0;
         rshift = 8;
         gshift = 16;
         bshift = 24;
      }
      else {
         ashift = 24;
         rshift = 16;
         gshift = 8;
         bshift = 0;
      }
      rgbmode = GL_TRUE;
   }
   else if (format==OSMESA_RGB) {
      indexBits = 0;
      redBits = CHAN_BITS;
      greenBits = CHAN_BITS;
      blueBits = CHAN_BITS;
      alphaBits = 0;
      bshift = 0;
      gshift = 8;
      rshift = 16;
      ashift = 24;
      rind = 0;
      gind = 1;
      bind = 2;
      rgbmode = GL_TRUE;
   }
   else if (format==OSMESA_BGR) {
      indexBits = 0;
      redBits = CHAN_BITS;
      greenBits = CHAN_BITS;
      blueBits = CHAN_BITS;
      alphaBits = 0;
      bshift = 0;
      gshift = 8;
      rshift = 16;
      ashift = 24;
      rind = 2;
      gind = 1;
      bind = 0;
      rgbmode = GL_TRUE;
   }
#if CHAN_TYPE == GL_UNSIGNED_BYTE
   else if (format==OSMESA_RGB_565) {
      indexBits = 0;
      redBits = 5;
      greenBits = 6;
      blueBits = 5;
      alphaBits = 0;
      rshift = 11;
      gshift = 5;
      bshift = 0;
      ashift = 0;
      rind = 0; // not used
      gind = 0;
      bind = 0;
      rgbmode = GL_TRUE;
   }
#endif
   else {
      return NULL;
   }


   osmesa = (OSMesaContext) CALLOC_STRUCT(osmesa_context);
   if (osmesa) {
      osmesa->gl_visual = _mesa_create_visual( rgbmode,
                                               GL_FALSE,    // double buffer
                                               GL_FALSE,    // stereo
                                               redBits,
                                               greenBits,
                                               blueBits,
                                               alphaBits,
                                               indexBits,
                                               depthBits,
                                               stencilBits,
                                               accumBits,
                                               accumBits,
                                               accumBits,
                                               alphaBits ? accumBits : 0,
                                               1            // num samples
                                               );
      if (!osmesa->gl_visual) {
         FREE(osmesa);
         return NULL;
      }

      // Initialize device driver function table 
      _mesa_init_driver_functions(&functions);
      // override with our functions
      functions.GetString = get_string;
      functions.UpdateState = osmesa_update_state;
      functions.GetBufferSize = get_buffer_size;
      functions.Clear = clear;

      if (!_mesa_initialize_context(&osmesa->mesa,
                                    osmesa->gl_visual,
                                    sharelist ? &sharelist->mesa
                                              : (GLcontext *) NULL,
                                    &functions, (void *) osmesa)) {
         _mesa_destroy_visual( osmesa->gl_visual );
         FREE(osmesa);
         return NULL;
      }

      _mesa_enable_sw_extensions(&(osmesa->mesa));
      _mesa_enable_1_3_extensions(&(osmesa->mesa));
      _mesa_enable_1_4_extensions(&(osmesa->mesa));
      _mesa_enable_1_5_extensions(&(osmesa->mesa));

      osmesa->gl_buffer = _mesa_create_framebuffer( osmesa->gl_visual,
                           (GLboolean) ( osmesa->gl_visual->depthBits > 0 ),
                           (GLboolean) ( osmesa->gl_visual->stencilBits > 0 ),
                           (GLboolean) ( osmesa->gl_visual->accumRedBits > 0 ),
                           GL_FALSE ); // s/w alpha

      if (!osmesa->gl_buffer) {
         _mesa_destroy_visual( osmesa->gl_visual );
         _mesa_free_context_data( &osmesa->mesa );
         FREE(osmesa);
         return NULL;
      }
      osmesa->format = format;
      osmesa->buffer = NULL;
      osmesa->width = 0;
      osmesa->height = 0;
      osmesa->userRowLength = 0;
      osmesa->rowlength = 0;
      osmesa->yup = GL_TRUE;
      osmesa->rshift = rshift;
      osmesa->gshift = gshift;
      osmesa->bshift = bshift;
      osmesa->ashift = ashift;
      osmesa->rInd = rind;
      osmesa->gInd = gind;
      osmesa->bInd = bind;
      osmesa->aInd = aind;

      // Initialize the software rasterizer and helper modules/
      {
	 GLcontext *ctx = &osmesa->mesa;

	 if (!_swrast_CreateContext( ctx ) ||
             !_ac_CreateContext( ctx ) ||
             !_tnl_CreateContext( ctx ) ||
             !_swsetup_CreateContext( ctx )) {
            _mesa_destroy_visual(osmesa->gl_visual);
            _mesa_free_context_data(ctx);
            _mesa_free(osmesa);
            return NULL;
         }
	
	 _swsetup_Wakeup( ctx );
         hook_in_driver_functions( ctx );
      }
   }
*/
   return osmesa;
}

//-------------------------------------------------------------------------
/*
 * Bind an OSMesaContext to an image buffer.  The image buffer is just a
 * block of memory which the client provides.  Its size must be at least
 * as large as width*height*sizeof(type).  Its address should be a multiple
 * of 4 if using RGBA mode.
 *
 * Image data is stored in the order of glDrawPixels:  row-major order
 * with the lower-left image pixel stored in the first array position
 * (ie. bottom-to-top).
 *
 * If the context's viewport hasn't been initialized yet, it will now be
 * initialized to (0,0,width,height).
 *
 * Input:  ctx - the rendering context
 *         buffer - the image buffer memory
 *         type - data type for pixel components
 *            Normally, only GL_UNSIGNED_BYTE and GL_UNSIGNED_SHORT_5_6_5
 *            are supported.  But if Mesa's been compiled with CHAN_BITS==16
 *            then type must be GL_UNSIGNED_SHORT.  And if Mesa's been build
 *            with CHAN_BITS==32 then type must be GL_FLOAT.
 *         width, height - size of image buffer in pixels, at least 1
 * Return:  GL_TRUE if success, GL_FALSE if error because of invalid ctx,
 *          invalid buffer address, invalid type, width<1, height<1,
 *          width>internal limit or height>internal limit.
 */
GLboolean MOSMesaMakeCurrent(MOSMesaContext ctx, void *buffer, GLenum type, GLsizei width, GLsizei height)
{

/**
   if (!ctx || !buffer ||
       width < 1 || height < 1 ||
       width > MAX_WIDTH || height > MAX_HEIGHT) {
      return GL_FALSE;
   }

   if (ctx->format == OSMESA_RGB_565) {
      if (type != GL_UNSIGNED_SHORT_5_6_5)
         return GL_FALSE;
   }
   else if (type != CHAN_TYPE) {
      return GL_FALSE;
   }

   osmesa_update_state( &ctx->mesa, 0 );
   _mesa_make_current( &ctx->mesa, ctx->gl_buffer );

   ctx->buffer = buffer;
   ctx->width = width;
   ctx->height = height;
   if (ctx->userRowLength)
      ctx->rowlength = ctx->userRowLength;
   else
      ctx->rowlength = width;

   compute_row_addresses( ctx );

   // init viewport 
   if (ctx->mesa.Viewport.Width == 0) {
      // initialize viewport and scissor box to buffer size 
      _mesa_Viewport( 0, 0, width, height );
      ctx->mesa.Scissor.Width = width;
      ctx->mesa.Scissor.Height = height;
   }
   else {
      // this will make ensure we recognize the new buffer size 
      _mesa_ResizeBuffersMESA();
   }

   // Added by Gerk Huisma: 
   _tnl_MakeCurrent( &ctx->mesa, ctx->mesa.DrawBuffer,
                     ctx->mesa.ReadBuffer );
**/
   return GL_TRUE;
}

//-------------------------------------------------------------------------
