/*	This file contains utilities for XOPs that create windows.
	It includes utilities for IGOR text windows (TU windows).
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

/*	GetActiveWindowRef()
	
	Returns an XOP_WINDOW_REF for the active window. An XOP_WINDOW_REF is a
	WindowPtr on Macintosh and an HWND on Windows. The returned value could be
	NULL and it could be a reference to a window that is not owned by the calling
	XOP.	
	
	Thread Safety: GetActiveWindowRef is not thread-safe.
*/
XOP_WINDOW_REF
GetActiveWindowRef(void)
{
	if (!CheckRunningInMainThread("GetActiveWindowRef"))
		return NULL;

	#ifdef MACIGOR
		return FrontNonFloatingWindow();
	#endif
	#ifdef WINIGOR
	{
		HWND hwndActive;
		
		/*	GetActiveWindows does not work when the window is an MDI child.
			The Windows documentation is characteristically unclear about this.
		*/
		hwndActive = (HWND)SendMessage(IgorClientHWND(), WM_MDIGETACTIVE, 0, 0);
		return hwndActive;
	}
	#endif
}

/*	IsXOPWindowActive(windowRef)

	Thread Safety: IsXOPWindowActive is not thread-safe.
*/
int
IsXOPWindowActive(XOP_WINDOW_REF windowRef)
{
	if (!CheckRunningInMainThread("IsXOPWindowActive"))
		return 0;

	if (windowRef == NULL)
		return 0;
	if (windowRef == GetActiveWindowRef())
		return 1;
	return 0;
}

/*	ShowXOPWindow(windowRef)

	Shows the window without activating it. Call this in response to the XOP_SHOW_WINDOW
	message from Igor.
	
	Thread Safety: ShowXOPWindow is not thread-safe.
*/
void
ShowXOPWindow(XOP_WINDOW_REF windowRef)
{
	if (!CheckRunningInMainThread("ShowXOPWindow"))
		return;

	#ifdef MACIGOR
		ShowWindow(windowRef);
	#endif
	#ifdef WINIGOR
		ShowWindow(windowRef, SW_SHOW);
	#endif
}

/*	HideXOPWindow(windowRef)

	Hides the window without sending it to the bottom of the desktop. Call this in response
	to the XOP_SHOW_WINDOW message from Igor.
	
	Thread Safety: HideXOPWindow is not thread-safe.
*/
void
HideXOPWindow(XOP_WINDOW_REF windowRef)
{
	if (!CheckRunningInMainThread("HideXOPWindow"))
		return;

	#ifdef MACIGOR
		HideWindow(windowRef);
	#endif
	#ifdef WINIGOR
		ShowWindow(windowRef, SW_HIDE);
	#endif
}

/*	ShowAndActivateXOPWindow(windowRef)

	This routine was added in XOP Toolkit release 3.1 but works with all supported
	versions of Igor.
	
	Thread Safety: ShowAndActivateXOPWindow is not thread-safe.
*/
void
ShowAndActivateXOPWindow(XOP_WINDOW_REF windowRef)
{
	if (!CheckRunningInMainThread("ShowAndActivateXOPWindow"))
		return;

	#ifdef MACIGOR
		ShowWindow(windowRef);
		SelectWindow(windowRef);
	#endif
	#ifdef WINIGOR
		ShowWindow(windowRef, SW_SHOW);
		SetWindowPos(windowRef, HWND_TOP, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
	#endif
}

/*	HideAndDeactivateXOPWindow(windowRef)

	This routine was added in XOP Toolkit release 3.1 but works with all supported
	versions of Igor.
	
	Thread Safety: HideAndDeactivateXOPWindow is not thread-safe.
*/
void
HideAndDeactivateXOPWindow(XOP_WINDOW_REF windowRef)
{
	if (!CheckRunningInMainThread("HideAndDeactivateXOPWindow"))
		return;

	#ifdef MACIGOR
		HideWindow(windowRef);
		SendBehind(windowRef, NULL);
	#endif
	#ifdef WINIGOR
	{
		HWND nextHWND;
		int isActive;
		
		nextHWND = GetWindow(windowRef, GW_HWNDNEXT);
		isActive = IsXOPWindowActive(windowRef);
		ShowWindow(windowRef, SW_HIDE);
		if (isActive) {
			// SetActiveWindow(nextHWND);	// Does not work when hwnd is MDI child, although it should, according to Windows documentation.
			SetWindowPos(nextHWND, HWND_TOP, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
		}
	}
	#endif
}

/*	SetXOPWindowTitle(windowRef,  title)

	This routine was added in XOP Toolkit release 3.1 but works with all supported
	versions of Igor.
	
	Thread Safety: SetXOPWindowTitle is not thread-safe.
*/
void
SetXOPWindowTitle(XOP_WINDOW_REF windowRef, const char* title)
{
	if (!CheckRunningInMainThread("SetXOPWindowTitle"))
		return;

	#ifdef MACIGOR
	{
		unsigned char pTitle[256];
		CopyCStringToPascal(title, pTitle);
		SetWTitle(windowRef, pTitle);
	}
	#endif
	#ifdef WINIGOR
		SetWindowText(windowRef, title);
	#endif
}

/*	GetXOPWindowPositionAndState(theWindow, r, winStatePtr)
	
	Returns the XOP window's position on the screen in pixels and its state.
	Used with SetXOPWindowPositionAndState to save and restore a window's position
	and state.
	
	Use this routine when you need to store a window position in a platform-dependent
	way, for example, in a preference file. Use GetXOPWindowIgorPositionAndState
	to store a window position in a platform-independent way, for example, in a
	/W=(left,top,right,bottom) flag.
	
	On Macintosh, the returned coordinates specify the location of the window's
	content region in global coordinates. Bit 0 of *winStatePtr is set if the
	window is visible and cleared if it is hidden. All other bits are set to 0.
	
	On Windows, the returned coordinates specify the the location of the entire
	window in its normal state relative the the top/left corner of the Igor MDI
	client window. Bit 0 of *winStatePtr is set if the window is visible and
	cleared if it is hidden. Bit 1 of *winStatePtr is set if the window is minimize
	and cleared if it is not minimized. All other bits are set to 0.
	
	On either platform, the returned rectangle is a Macintosh rectangle.

	This routine was added in XOP Toolkit release 3.1 but works with all supported
	versions of Igor.
	
	Thread Safety: GetXOPWindowPositionAndState is not thread-safe.
*/
void
GetXOPWindowPositionAndState(XOP_WINDOW_REF windowRef, Rect* r, int* winStatePtr)
{
	if (!CheckRunningInMainThread("GetXOPWindowPositionAndState")) {
		MemClear(r, sizeof(Rect));
		*winStatePtr = 0;
		return;
	}

	#ifdef MACIGOR
	{
		GrafPtr thePort;
		GrafPtr savePort;

		thePort = GetWindowPort(windowRef);
		*winStatePtr = IsWindowVisible(windowRef) ? 1:0;
		GetPortBounds(thePort, r);				// Content rectangle in local coordinates.
		GetPort(&savePort);
		SetPort(thePort);
		LocalToGlobal((Point *)&r->top);
		LocalToGlobal((Point *)&r->bottom);
		SetPort(savePort);
	}
	#endif
	
	#ifdef WINIGOR
	{
		WINDOWPLACEMENT wp;
		
		*winStatePtr = IsWindowVisible(windowRef) ? 1:0;
		if (IsIconic(windowRef))
			*winStatePtr |= 2;
		wp.length = sizeof(wp);
		GetWindowPlacement(windowRef, &wp);
		WinRectToMacRect(&wp.rcNormalPosition, r);
	}
	#endif
}

/*	SetXOPWindowPositionAndState(theWindow, r, winState)
	
	Moves the XOP window to the position indicated by r and sets its state.
	Used with GetXOPWindowPositionAndState to save and restore a window's position
	and state.
	
	Use this routine when you need to restore a window position in a platform-dependent
	way, for example, in a preference file. Use SetXOPWindowIgorPositionAndState
	to restore a window position in a platform-independent way, for example, in a
	/W=(left,top,right,bottom) flag.
	
	See GetXOPWindowPositionAndState for a discussion of the units of the rectangle
	and the meaning of the winState parameter.
	
	This routine makes an effort to prevent the window from becoming inaccessible
	because it is off-screen.

	This routine was added in XOP Toolkit release 3.1 but works with all supported
	versions of Igor.
	
	Thread Safety: SetXOPWindowPositionAndState is not thread-safe.
*/
void
SetXOPWindowPositionAndState(XOP_WINDOW_REF windowRef, Rect* r, int winState)
{
	if (!CheckRunningInMainThread("SetXOPWindowPositionAndState"))
		return;

	#ifdef MACIGOR
	{
		int left, right, top, bottom, width, height;
		RgnHandle grayRgn;
		Rect frameRect;

		if ((winState & 1) == 0)
			HideWindow(windowRef);

		left = r->left;
		right = r->right;
		top = r->top;
		bottom = r->bottom;
		width = right - left;
		height = bottom - top;
		
		SetRect(&frameRect, left, top-20, right, top);		// Roughly the window frame rect.
		InsetRect(&frameRect, 5, 5);						// Make sure at least 5 pixels of frame are on the desktop.
		grayRgn = GetGrayRgn();
		if (!RectInRgn(&frameRect, grayRgn)) {				// Window frame is completely out of desktop?
			Rect bounds;
			GetRegionBounds(grayRgn, &bounds);
			left = bounds.left + 3;
			top = bounds.top + 23;
			right = left + width;
			bottom = top + height;
		}
		if (width < 10)
			width = 10;
		if (height < 10)
			height = 10;

		MoveWindow(windowRef, left, top, 0);
		SizeWindow(windowRef, width, height, -1);

		if ((winState & 1) != 0)
			ShowWindow(windowRef);
	}
	#endif
	
	#ifdef WINIGOR
	{
		WINDOWPLACEMENT wp;
		RECT childRECT, mdiClientRECT;

		if ((winState & 1) == 0)
			ShowWindow(windowRef, SW_HIDE);
		
		wp.length = sizeof(wp);
		GetWindowPlacement(windowRef, &wp);
		wp.flags = 0;
		if (winState & 2) {						// Want to minimize?
			wp.showCmd = SW_MINIMIZE;
		}
		else {
			if ((winState & 1) != 0)
				wp.showCmd = SW_SHOWNORMAL;
			else
				wp.showCmd = SW_HIDE;
		}
		
		MacRectToWinRect(r, &childRECT);
		GetClientRect(IgorClientHWND(), &mdiClientRECT);
		if (childRECT.top < 0)
			OffsetRect(&childRECT, 0, -childRECT.top);
		if (childRECT.top > mdiClientRECT.bottom)
			OffsetRect(&childRECT, 0, -(childRECT.top - mdiClientRECT.bottom + 10));
		if (childRECT.right < 0)
			OffsetRect(&childRECT, -childRECT.right + 10, 0);
		if (childRECT.left > mdiClientRECT.right)
			OffsetRect(&childRECT, -(childRECT.left - mdiClientRECT.right + 10), 0);
		wp.rcNormalPosition = childRECT;
		
		SetWindowPlacement(windowRef, &wp);
	}
	#endif
}

/*	TransformWindowCoordinates(mode, coords)

	Transforms window coordinates from screen pixels into Igor coordinates or
	from Igor coordinates into screen pixels. This routine is intended for use
	in command line operations that set a window position, for example, for
	an operation that supports a /W=(left,top,right,bottom) flag. We want
	a given command containing a /W flag to produce approximately the same result
	on Macintosh and on Windows. This is complicated because of differences in
	the way each platform represents the position of windows.

	Igor coordinates are a special kind of coordinates that were designed to solve
	this problem. Igor coordinates are in units of points, not pixels. On Macintosh,
	Igor coordinates are the same as global coordinates - points relative to the top/left
	corner of the main screen. On Windows, Igor coordinates are points relative to a spot
	20 points above the top/left corner of the MDI client area. As a result of this
	definition, the vertical coordinate 20 corresponds to the point just below the main
	menu bar on both platforms.

	The use of Igor coordinates in commands means that you can transport files and
	commands from one platform to the other and get reasonable results. However,
	the results may not be exactly what you expect. The reason for this is that
	Igor positions the "content" portion of a window. The content portion is the
	portion of the window exclusive of the frame and title bar (border and caption
	in Windows terminology). Because the size of window borders and captions is
	variable on Windows, when you open a Macintosh experiment on Windows or vice versa,
	the window positions might be slightly different from what you would expect. 
	
	We keep coordinates in floating point because, to accurately reposition a
	window, we need to use fractional points in /W=(left,top,right,bottom) flags.
	
	mode is
		0:		Transform from screen pixels into Igor coordinates.
		1:		Transform from Igor coordinates into screen pixels.
		
	For TransformWindowCoordinates, screen pixels are in global coordinates on
	Macintosh (relative to the top/left corner of the main screen) and are in
	MDI-client coordinates (relative to the top/left corner of the MDI client
	window, not the MDI frame) on Windows.  
		
	coords is an array of window coordinates. It is both an input and an output.
	The coordinates specify the location of the window's content area only. That
	is, it excludes the title bar and the frame.
	
	coords[0] is the location of the left edge of the window content area.
	coords[1] is the location of the top edge of the window content area.
	coords[2] is the location of the right edge of the window content area.
	coords[3] is the location of the bottom edge of the window content area.
	
	On Macintosh, screen coordinates and Igor coordinates are identical. Thus, this
	routine is a NOP on Macintosh.

	This routine was added in XOP Toolkit release 3.1 but works with all supported
	versions of Igor.
	
	Thread Safety: TransformWindowCoordinates is not thread-safe.
*/
void
TransformWindowCoordinates(int mode, double coords[4])
{
	if (!CheckRunningInMainThread("TransformWindowCoordinates"))
		return;

	#ifdef MACIGOR			// [
		#pragma unused(mode)
		#pragma unused(coords)
	#endif					// ]

	#ifdef WINIGOR			// [
	{
		HDC screenDC;
		int hPixelsPerInch, vPixelsPerInch;
		int vOffsetInPoints;
	
		screenDC = GetDC(NULL);
		hPixelsPerInch = GetDeviceCaps(screenDC,LOGPIXELSX);
		vPixelsPerInch = GetDeviceCaps(screenDC,LOGPIXELSY);
		ReleaseDC(NULL,screenDC);
	
		vOffsetInPoints = 20;
		
		if (mode == 0) {
			// Transform from screen pixels into IGOR coordinates.
			coords[0] = 72.0 * (coords[0] / hPixelsPerInch);
			coords[1] = vOffsetInPoints + 72.0 * (coords[1] / vPixelsPerInch);
			coords[2] = 72.0 * (coords[2] / hPixelsPerInch);
			coords[3] = vOffsetInPoints + 72.0 * (coords[3] / vPixelsPerInch);
		}
		
		if (mode == 1) {
			// Transform from IGOR coordinates into screen pixels.
			coords[0] = (coords[0] / 72.0) * hPixelsPerInch;
			coords[1] = ((coords[1] - vOffsetInPoints) / 72.0) * vPixelsPerInch;
			coords[2] = (coords[2] / 72.0) * hPixelsPerInch;
			coords[3] = ((coords[3] - vOffsetInPoints) / 72.0) * vPixelsPerInch;
		}
	}
	#endif					// ]
}

/*	GetXOPWindowIgorPositionAndState(theWindow, coords, winStatePtr)
	
	Returns the XOP window's position on the screen in Igor coordinates and its state.
	Used with SetXOPWindowIgorPositionAndState to save and restore a window's position
	and state.
	
	Use this routine when you need to store a window position in a platform-independent
	way, for example, in a /W=(left,top,right,bottom) flag. Use GetXOPWindowPositionAndState
	to store a window position in a platform-dependent way, for example, in a preference file.
	
	See TransformWindowCoordinates for a discussion of Igor coordinates. On
	both Macintosh and Windows, the returned coordinates specify the location
	of the window's content region, not the outside edges of the window. On Windows,
	the returned coordinates specify the the location of the window in its normal state
	even if the window is minmized or maximized.
	
	On Macintosh, bit 0 of *winStatePtr is set if the window is visible and cleared
	if it is hidden. All other bits are set to 0.
	
	On Windows, bit 0 of *winStatePtr is set if the window is visible and
	cleared if it is hidden. Bit 1 of *winStatePtr is set if the window is minimize
	and cleared if it is not minimized. All other bits are set to 0.

	This routine was added in XOP Toolkit release 3.1 but works with all supported
	versions of Igor.
	
	Thread Safety: GetXOPWindowIgorPositionAndState is not thread-safe.
*/
void
GetXOPWindowIgorPositionAndState(XOP_WINDOW_REF windowRef, double coords[4], int* winStatePtr)
{
	if (!CheckRunningInMainThread("GetXOPWindowIgorPositionAndState")) {
		MemClear(coords, 4*sizeof(double));
		*winStatePtr = 0;
		return;
	}

	#ifdef MACIGOR
	{
		Rect r;
		GrafPtr thePort;
		GrafPtr savePort;

		thePort = GetWindowPort(windowRef);
		*winStatePtr = IsWindowVisible(windowRef) ? 1:0;
		GetPortBounds(thePort, &r);				// Content rectangle in local coordinates.
		GetPort(&savePort);
		SetPort(thePort);
		LocalToGlobal((Point *)&r.top);
		LocalToGlobal((Point *)&r.bottom);
		SetPort(savePort);
		coords[0] = r.left;
		coords[1] = r.top;
		coords[2] = r.right;
		coords[3] = r.bottom;
		TransformWindowCoordinates(0, coords);
	}
	#endif
	
	#ifdef WINIGOR
	{
		WINDOWPLACEMENT wp;
		RECT wr;
		int frameHeight, frameWidth, captionHeight;
		
		*winStatePtr = IsWindowVisible(windowRef) ? 1:0;
		if (IsIconic(windowRef))
			*winStatePtr |= 2;
		wp.length = sizeof(wp);
		GetWindowPlacement(windowRef, &wp);
		wr = wp.rcNormalPosition;
		
		// wr contains the outer window coordinates. We transform this into the "content area".
		frameHeight = GetSystemMetrics(SM_CYFRAME);
		frameWidth = GetSystemMetrics(SM_CXFRAME);
		captionHeight = GetSystemMetrics(SM_CYCAPTION);
		wr.left += frameWidth;
		wr.right -= frameWidth;
		wr.top += frameHeight + captionHeight;
		wr.bottom -= frameHeight;
	
		// Now convert to floating point and to Igor coordinates.
		coords[0] = wr.left;
		coords[1] = wr.top;
		coords[2] = wr.right;
		coords[3] = wr.bottom;
		TransformWindowCoordinates(0, coords);
	}
	#endif
}

/*	SetXOPWindowIgorPositionAndState(theWindow, coords, winState)
	
	Moves the XOP window to the position indicated by coords and sets its state.
	Used with GetXOPWindowIgorPositionAndState to save and restore a window's position
	and state.
	
	Use this routine when you need to restore a window position in a platform-independent
	way, for example, in a /W=(left,top,right,bottom) flag. Use SetXOPWindowPositionAndState
	to restore a window position in a platform-dependent way, for example, in a preference file.
	
	See GetXOPWindowIgorPositionAndState for a discussion of the units of coords
	and the meaning of the winState parameter.
	
	This routine makes an effort to prevent the window from becoming inaccessible
	because it is off-screen.

	This routine was added in XOP Toolkit release 3.1 but works with all supported
	versions of Igor.
	
	Thread Safety: SetXOPWindowIgorPositionAndState is not thread-safe.
*/
void
SetXOPWindowIgorPositionAndState(XOP_WINDOW_REF windowRef, double coords[4], int winState)
{
	if (!CheckRunningInMainThread("SetXOPWindowIgorPositionAndState"))
		return;

	#ifdef MACIGOR
	{
		int left, right, top, bottom, width, height;
		RgnHandle grayRgn;
		Rect frameRect;
		double coords2[4];
		
		memcpy(coords2, coords, sizeof(coords2));			// To avoid modifying the input parameter.
		TransformWindowCoordinates(1, coords2);				// Transform from Igor coordinates into screen coordinates.

		if ((winState & 1) == 0)
			HideWindow(windowRef);

		left = coords2[0] + .5;								// + .5 for rounding.
		right = coords2[2] + .5;
		top = coords2[1] + .5;
		bottom = coords2[3] + .5;
		width = right - left;
		height = bottom - top;
		
		SetRect(&frameRect, left, top-20, right, top);		// Roughly the window frame rect.
		InsetRect(&frameRect, 5, 5);						// Make sure at least 5 pixels of frame are on the desktop.
		grayRgn = GetGrayRgn();
		if (!RectInRgn(&frameRect, grayRgn)) {				// Window frame is completely out of desktop?
			Rect bounds;
			GetRegionBounds(grayRgn, &bounds);
			left = bounds.left + 3;
			top = bounds.top + 23;
			right = left + width;
			bottom = top + height;
		}
		if (width < 10)
			width = 10;
		if (height < 10)
			height = 10;

		MoveWindow(windowRef, left, top, 0);
		SizeWindow(windowRef, width, height, -1);

		if ((winState & 1) != 0)
			ShowWindow(windowRef);
	}
	#endif
	
	#ifdef WINIGOR
	{
		WINDOWPLACEMENT wp;
		RECT childRECT, mdiClientRECT;
		double coords2[4];
		int frameHeight, frameWidth, captionHeight;
		
		memcpy(coords2, coords, sizeof(coords2));			// To avoid modifying the input parameter.
		TransformWindowCoordinates(1, coords2);				// Transform from Igor coordinates into screen coordinates.

		if ((winState & 1) == 0)
			ShowWindow(windowRef, SW_HIDE);
		
		wp.length = sizeof(wp);
		GetWindowPlacement(windowRef, &wp);
		wp.flags = 0;
		if (winState & 2) {						// Want to minimize?
			wp.showCmd = SW_MINIMIZE;
		}
		else {
			if ((winState & 1) != 0)
				wp.showCmd = SW_SHOWNORMAL;
			else
				wp.showCmd = SW_HIDE;
		}
		
		// HR, 091001: Added casts to prevent VC warnings.
		childRECT.left = (LONG)(coords2[0] + .5);					// + .5 for rounding.
		childRECT.right = (LONG)(coords2[2] + .5);
		childRECT.top = (LONG)(coords2[1] + .5);
		childRECT.bottom = (LONG)(coords2[3] + .5);
		
		// childRECT contains the "content area" coordinates. We transform this into the outer window coordinates.
		frameHeight = GetSystemMetrics(SM_CYFRAME);
		frameWidth = GetSystemMetrics(SM_CXFRAME);
		captionHeight = GetSystemMetrics(SM_CYCAPTION);
		childRECT.left -= frameWidth;
		childRECT.right += frameWidth;
		childRECT.top -= frameHeight + captionHeight;
		childRECT.bottom += frameHeight;

		GetClientRect(IgorClientHWND(), &mdiClientRECT);
		if (childRECT.top < 0)
			OffsetRect(&childRECT, 0, -childRECT.top);
		if (childRECT.top > mdiClientRECT.bottom)
			OffsetRect(&childRECT, 0, -(childRECT.top - mdiClientRECT.bottom + 10));
		if (childRECT.right < 0)
			OffsetRect(&childRECT, -childRECT.right + 10, 0);
		if (childRECT.left > mdiClientRECT.right)
			OffsetRect(&childRECT, -(childRECT.left - mdiClientRECT.right + 10), 0);
		wp.rcNormalPosition = childRECT;
		
		SetWindowPlacement(windowRef, &wp);
	}
	#endif
}

/*	TellIgorWindowStatus(windowRef, status, options)

	Call TellIgorWindowStatus when you hide or show your window, when it is activated and deactivated
	and when it is about to be killed and was killed.
	
	TellIgorWindowStatus allows your window to be included various Igor features such
	as the Recent Windows menu item and the Show Recently Hidden Windows menu item.
	It is not required that you call this routine but is recommended. If you do not call
	this routine then your windows will not participate in these features.
	
	Here are the calls to make at the appropriate times:
		TellIgorWindowStatus(windowRef, WINDOW_STATUS_DID_HIDE, 0);
		TellIgorWindowStatus(windowRef, WINDOW_STATUS_DID_SHOW, 0);
		TellIgorWindowStatus(windowRef, WINDOW_STATUS_ACTIVATED, 0);
		TellIgorWindowStatus(windowRef, WINDOW_STATUS_DEACTIVATED, 0);
		TellIgorWindowStatus(windowRef, WINDOW_STATUS_ABOUT_TO_KILL, 0);
		TellIgorWindowStatus(windowRef, WINDOW_STATUS_DID_KILL, 0);			// Added in Igor Pro 6.04B01.
	
	NOTE: If you call this routine with any status, you must also call it when your window
	is about to be killed so Igor can remove it from various internal lists. Failure to do
	so may cause a crash.
	
	Returns 0 if OK or IGOR_OBSOLETE or NOT_IN_THREADSAFE.
	
	Thread Safety: TellIgorWindowStatus is not thread-safe.
*/
int
TellIgorWindowStatus(XOP_WINDOW_REF windowRef, int status, int options)
{
	if (!CheckRunningInMainThread("TellIgorWindowStatus"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(TELL_IGOR_WINDOW_STATUS, windowRef, (void*)status, (void*)options);
}


/*	*** Text Windows for XOPs *** */

// TUNew is obsolete and was removed from XOP Toolkit 6. Use TUNew2.

/*	TUNew2(winTitle, winRectPtr, TUPtr, windowRefPtr)

	TUNew2 creates a new text window and TU document.
	
	winTitle points to the title for the new window.
	
	winRectPtr points to a Macintosh Rect which specifies the size and location of
	the content region of the window in units of pixels.
	
	On Macintosh, this is in global coordinates.  Use a top coordinate of 40 to position
	the window right below the menu bar.
	
	On Windows, it is in client window coordinates of the Igor MDI frame window. Use a
	top coordinate of 22 to position the window right below the menu bar.
	
	It returns via TUPtr a handle to the TU document and returns via windowRefPtr
	a pointer to a WindowPtr (Mac) or HWND (Windows) for the newly created window.
	
	In the event of an error, it returns non-zero as the function result and NULL via
	TUPtr and windowRefPtr.

	TUNew2 uses a default font and font size. The resulting text document is like
	an Igor plain text notebook.
	
	The window is initially hidden. Call ShowAndActivateXOPWindow to show it.
	
	Thread Safety: TUNew2 is not thread-safe.
*/
int
TUNew2(const char* winTitle, const Rect* winRectPtr, Handle* TUPtr, XOP_WINDOW_REF* windowRefPtr)
{
	int err;

	*TUPtr = NULL;
	*windowRefPtr = NULL;

	if (!CheckRunningInMainThread("TUNew2"))
		return NOT_IN_THREADSAFE;

	err = (int)CallBack4(TUNEW2, (void*)winTitle, (void*)winRectPtr, TUPtr, windowRefPtr);
	return err;
}
	
/*	TUDispose(TU)

	TUDispose disposes of the TE record, the window, and the TU record for the specified
	TU window. This should be called for any TU windows when the XOP is about to be disposed.
	
	Thread Safety: TUDispose is not thread-safe.
*/
void	
TUDispose(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TUDispose"))
		return;

	CallBack1(TUDISPOSE, TU);
}
	
/*	TUDisplaySelection(TU)

	Tries to get the selected text in view as best as it can by scrolling.
	
	The rules for vertical scrolling used by TUSelToView are:
		if selStart and selEnd are in view, do nothing
		if selected text won't fit in window vertically
			bring selStart line to top of window
		if selected text will fit vertically
			if selStart line is above
				bring selStart line to top
			if selEnd is below
				bring selEnd line to bottom
	
	If the selected text is multiline, it won't scroll horizontally to get the right edge in
	view.  This gives best intuitive results.
	
	Thread Safety: TUDisplaySelection is not thread-safe.
*/
void
TUDisplaySelection(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TUDisplaySelection"))
		return;

	CallBack1(TUDISPLAYSELECTION, TU);
}

/*	TUGrow(TU, size)
	
	TUGrow() adjust the window size.

	Size is the size of the window packed into a 32 bit integer.
	
	The vertical size is in the high word and the horizontal size is in the low word.
	
	However, if size = 0 then it does a zoom rather than a grow.
	
	Also, if size == -1, then TUGrow does not resize the window but merely adjusts
	for a change in size that has already been done. For example, TUGrow moves
	the scroll bars to the new edges of the window.
	
	Thread Safety: TUGrow is not thread-safe.
*/
void
TUGrow(TUStuffHandle TU, int size)
{	
	if (!CheckRunningInMainThread("TUGrow"))
		return;

	CallBack2(TUGROW, TU, (void*)size);
}

/*	TUDrawWindow(TU)

	Draws the window containing the text referred to by TU.
	
	Thread Safety: TUDrawWindow is not thread-safe.
*/
void
TUDrawWindow(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUDrawWindow"))
		return;

	CallBack1(TUDRAWWINDOW, TU);
}

/*	TUUpdate(TU)

	Updates the window containing the text referred to by TU if the updateRgn of the window
	is not empty.
	
	Thread Safety: TUUpdate is not thread-safe.
*/
void
TUUpdate(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUUpdate"))
		return;

	CallBack1(TUUPDATE, TU);
}

/*	TUFind(TU, code)

	code:	1 = normal find
			2 = find same
			3 = find selection
	
	Thread Safety: TUFind is not thread-safe.
*/
void
TUFind(TUStuffHandle TU, int code)
{	
	if (!CheckRunningInMainThread("TUFind"))
		return;

	CallBack2(TUFIND, TU, (void*)code);
}

/*	TUReplace(TU)
	
	Thread Safety: TUReplace is not thread-safe.
*/
void
TUReplace(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUReplace"))
		return;

	CallBack1(TUREPLACE, TU);
}

/*	TUIndentLeft(TU)
	
	Thread Safety: TUIndentLeft is not thread-safe.
*/
void
TUIndentLeft(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUIndentLeft"))
		return;

	CallBack1(TUINDENTLEFT, TU);
}

/*	TUIndentRight(TU)
	
	Thread Safety: TUIndentRight is not thread-safe.
*/
void
TUIndentRight(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUIndentRight"))
		return;

	CallBack1(TUINDENTRIGHT, TU);
}

/*	TUClick(TU, eventPtr)

	Services click referred to by eventPtr.
	
	Thread Safety: TUClick is not thread-safe.
*/
void
TUClick(TUStuffHandle TU, EventRecord* eventPtr)
{	
	if (!CheckRunningInMainThread("TUClick"))
		return;

	CallBack2(TUCLICK, TU, eventPtr);
}

/*	TUActivate(TU)
	
	Thread Safety: TUActivate is not thread-safe.
*/
void
TUActivate(TUStuffHandle TU, int flag)
{	
	if (!CheckRunningInMainThread("TUActivate"))
		return;

	CallBack2(TUACTIVATE, TU, (void*)flag);
}

/*	TUIdle(TU)
	
	Thread Safety: TUIdle is not thread-safe.
*/
void
TUIdle(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUIdle"))
		return;

	CallBack1(TUIDLE, TU);
}

/*	TUNull(TU)
	
	Thread Safety: TUNull is not thread-safe.
*/
void
TUNull(TUStuffHandle TU, EventRecord* eventPtr)
{	
	if (!CheckRunningInMainThread("TUNull"))
		return;

	CallBack2(TUNULL, TU, eventPtr);
}

/*	TUCopy(TU)
	
	Thread Safety: TUCopy is not thread-safe.
*/
void
TUCopy(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUCopy"))
		return;

	CallBack1(TUCOPY, TU);
}

/*	TUCut(TU)
	
	Thread Safety: TUCut is not thread-safe.
*/
void
TUCut(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUCut"))
		return;

	CallBack1(TUCUT, TU);
}

/*	TUPaste(TU)
	
	Thread Safety: TUPaste is not thread-safe.
*/
void
TUPaste(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUPaste"))
		return;

	CallBack1(TUPASTE, TU);
}

/*	TUClear(TU)
	
	Thread Safety: TUClear is not thread-safe.
*/
void
TUClear(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUClear"))
		return;

	CallBack1(TUCLEAR, TU);
}

/*	TUKey(TU)
	
	Thread Safety: TUKey is not thread-safe.
*/
void
TUKey(TUStuffHandle TU, EventRecord* eventPtr)
{	
	if (!CheckRunningInMainThread("TUKey"))
		return;

	CallBack3(TUKEY, TU, (void*)(int)eventPtr->modifiers, (void*)(eventPtr->message));
}

/*	TUInsert(TU)
	
	Thread Safety: TUInsert is not thread-safe.
*/
void
TUInsert(TUStuffHandle TU, const char *dataPtr, int dataLen)
{	
	if (!CheckRunningInMainThread("TUInsert"))
		return;

	CallBack3(TUINSERT, TU, (void*)dataPtr, (void*)dataLen);
}

/*	TUDelete(TU)
	
	Thread Safety: TUDelete is not thread-safe.
*/
void
TUDelete(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUDelete"))
		return;

	CallBack1(TUDELETE, TU);
}

// TUSetSelect is obsolete and was removed from XOP Toolkit 6. Use TUSetSelLocs instead.

/*	TUSelectAll(TU)
	
	Thread Safety: TUSelectAll is not thread-safe.
*/
void
TUSelectAll(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUSelectAll"))
		return;

	CallBack1(TUSELECTALL, TU);
}

/*	TUUndo(TU)
	
	Thread Safety: TUUndo is not thread-safe.
*/
void
TUUndo(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUUndo"))
		return;

	CallBack1(TUUNDO, TU);
}

/*	TUPrint(TU)
	
	Thread Safety: TUPrint is not thread-safe.
*/
void
TUPrint(TUStuffHandle TU)
{	
	if (!CheckRunningInMainThread("TUPrint"))
		return;

	CallBack1(TUPRINT, TU);
}

/*	TUFixEditMenu(TU)

	Sets items in edit menu properly according to the state of the TU document.
	
	Thread Safety: TUFixEditMenu is not thread-safe.
*/
void
TUFixEditMenu(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TUFixEditMenu"))
		return;

	CallBack1(TUFIXEDITMENU, TU);
}

/*	TUFixFileMenu(TU)

	Sets items in file menu properly according to the state of the TU document.
	
	Thread Safety: TUFixFileMenu is not thread-safe.
*/
void
TUFixFileMenu(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TUFixFileMenu"))
		return;

	CallBack1(TUFIXFILEMENU, TU);
}

// TUGetText is obsolete and was removed from XOP Toolkit 6. Use TUFetchParagraphText instead.

// TUFetchText is obsolete and was removed from XOP Toolkit 6. Use TUFetchParagraphText instead.

// TULength is obsolete and was removed from XOP Toolkit 6. Use TUGetDocInfo instead.

/*	TULines(TU)

	TULines returns the number of lines in the specified document.
	
	Thread Safety: TULines is not thread-safe.
*/
int
TULines(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TULines"))
		return 0;

	return (int)CallBack1(TULINES, TU);
}

// TUSelStart is obsolete and was removed from XOP Toolkit 6. Use TUGetSelLocs instead.

// TUSelEnd is obsolete and was removed from XOP Toolkit 6. Use TUGetSelLocs instead.

// TUSelectionLength is obsolete and was removed from XOP Toolkit 6. Use TUGetSelLocs instead.

// TUInsertFile is obsolete and was removed from XOP Toolkit 6. There is no direct replacement.

// TUWriteFile is obsolete and was removed from XOP Toolkit 6. There is no direct replacement.

/*	TUSFInsertFile(TU, prompt, fileTypes, numTypes)

	Gets file from user using standard file package open dialog.

	prompt is a prompt to appear in standard save dialog.
	fileTypes is a pointer to an OSType array of file types.
	numTypes is the number of file types in the array.
	
	Inserts text from the file at the insertion point of the specified document.
	Returns error code from insertion.
	
	Thread Safety: TUSFInsertFile is not thread-safe.
*/
int
TUSFInsertFile(TUStuffHandle TU, const char *prompt, OSType fileTypes[], int numTypes)
{	
	if (!CheckRunningInMainThread("TUSFInsertFile"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(TUSFINSERTFILE, TU, (void*)prompt, fileTypes, (void*)numTypes);
}

/*	TUSFWriteFile(TU, prompt, fileType, allFlag)

	Gets file name from user using standard file package save dialog.
	
	Writes text from the specified document to file. Replaces file if it already exists.
	
	prompt is a prompt to appear in standard save dialog.
	fileType is type of file to be written (e.g. 'TEXT').
	allFlag = 0 means write only selected text. Otherwise, write all text.
	
	Returns error code from write.
	
	Thread Safety: TUSFWriteFile is not thread-safe.
*/
int
TUSFWriteFile(TUStuffHandle TU, const char *prompt, OSType fileType, int allFlag)
{	
	if (!CheckRunningInMainThread("TUSFWriteFile"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(TUSFWRITEFILE, TU, (void*)prompt, (void*)fileType, (void*)allFlag);
}

/*	TUPageSetupDialog(TU)
	
	Thread Safety: TUPageSetupDialog is not thread-safe.
*/
void
TUPageSetupDialog(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TUPageSetupDialog"))
		return;

	CallBack1(TUPAGESETUPDIALOG, TU);
}

/*	TUGetDocInfo(TU, dip)
	
	Returns information about the text utility document via the structure pointed to by dip.
	You MUST execute
		dip->version = TUDOCINFO_VERSION
	before calling TUGetDocInfo so that Igor knows which version of the structure
	your XOP is using.
	
	Returns 0 if OK, -1 if unsupported version of the structure or an Igor error
	code if the version of Igor that is running does not support this callback.
	
	Thread Safety: TUGetDocInfo is not thread-safe.
*/
int
TUGetDocInfo(TUStuffHandle TU, TUDocInfoPtr dip)
{
	if (!CheckRunningInMainThread("TUGetDocInfo"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack2(TUGETDOCINFO, TU, dip);
}

/*	TUGetSelLocs(TU, startLocPtr, endLocPtr)

	Sets *startLocPtr and *endLocPtr to describe the selected text in the document.

	Returns 0 if OK, an Igor error code if the version of Igor that you are running
	with does not support this callback.
	
	Thread Safety: TUGetSelLocs is not thread-safe.
*/
int
TUGetSelLocs(TUStuffHandle TU, TULocPtr startLocPtr, TULocPtr endLocPtr)
{
	if (!CheckRunningInMainThread("TUGetSelLocs"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(TUGETSELLOCS, TU, startLocPtr, endLocPtr);
}

/*	TUSetSelLocs(TU, startLocPtr, endLocPtr, flags)

	If startLocPtr is not NULL, sets the selection in the text utility document
	based on startLocPtr and endLocPtr which must be valid.
	
	If flags is 1, displays the selection if it is out of view.
	Other bits in flags may be used for other purposes in the future.

	Returns 0 if OK, an Igor error code if the version of Igor that you are running
	with does not support this callback. Also returns an error if the start or
	end locations are out of bounds or if the start location is after the end location.
	
	Thread Safety: TUSetSelLocs is not thread-safe.
*/
int
TUSetSelLocs(TUStuffHandle TU, TULocPtr startLocPtr, TULocPtr endLocPtr, int flags)
{
	if (!CheckRunningInMainThread("TUSetSelLocs"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(TUSETSELLOCS, TU, startLocPtr, endLocPtr, (void*)flags);
}

/*	TUFetchParagraphText(TU, paragraph, textPtrPtr, lengthPtr)
	
	If textPtrPtr is not NULL, returns via textPtrPtr the text in the specified paragraph.
	
	Sets *lengthPtr to the number of characters in the paragraph whether textPtrPtr is NULL or not.
	
	paragraph is assumed to be a valid paragraph number.
	
	textPtrPtr is a pointer to your char* variable.
	Igor allocates a pointer, using NewPtr, and sets *textPtrPtr to point to the allocated memory.
	The returned pointer belongs to you. Dispose it using DisposePtr when you no longer need it.
	
	Example:
		char* p;
		int paragraph, length;
		int result;
		
		paragraph = 0;
		if (result = TUFetchParagraphText(TU, paragraph, &p, &length))
			return result;
		
		<Deal with the text>
		
		DisposePtr(p);
	
	Note that the text pointed to by p is NOT null terminated and therefore
	is not a C string. You can turn it into a C string as follows:
		SetPtrSize(p, length+1);
		if (result = MemError()) {
			DisposePtr(p);
			return result
		}
		p[length] = 0;

	Returns 0 if OK, an Igor error code if an error occurs fetching the text
	or the version of Igor that you are running with does not support this callback.
	Also returns an error if the paragraph is out of bounds.
	
	Thread Safety: TUFetchParagraphText is not thread-safe.
*/
int
TUFetchParagraphText(TUStuffHandle TU, int paragraph, Ptr *textPtrPtr, int *lengthPtr)
{
	if (!CheckRunningInMainThread("TUFetchParagraphText"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(TUFETCHPARAGRAPHTEXT, TU, (void*)paragraph, textPtrPtr, lengthPtr);
}

/*	TUFetchSelectedText(TU, textHandlePtr, reservedForFuture, flags)
	
	Returns via textHandlePtr the selected text in the text utility document.
	
	textHandlePtr is a pointer to your Handle variable.
	reservedForFuture should be NULL for now.
	flags should be 0 for now.
	
	Example:
		Handle h;
		int result;
		
		if (result = TUFetchSelectedText(TU, &h, NULL, 0))
			return result;
		
		<Deal with the text>
		
		DisposeHandle(h);
	
	Note that the text in the handle h is NOT null terminated and therefore
	is not a C string. You can turn it into a C string as follows:
		length = GetHandleSize(h);
		SetHandleSize(h, length+1);
		if (result = MemError()) {
			DisposeHandle(h);
			return result
		}
		*h[length] = 0;

	Returns 0 if OK, an Igor error code if an error occurs fetching the text
	or the version of Igor that you are running with does not support this callback.
	
	Thread Safety: TUFetchSelectedText is not thread-safe.
*/
int
TUFetchSelectedText(TUStuffHandle TU, Handle* textHandlePtr, void* reservedForFuture, int flags)
{
	if (!CheckRunningInMainThread("TUFetchSelectedText"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(TUFETCHSELECTEDTEXT, TU, textHandlePtr,  reservedForFuture, (void*)flags);
}

/*	TUFetchText2(TU, startLocPtr, endLocPtr, textHandlePtr, reservedForFuture, flags)
	
	Returns via textHandlePtr the text in the text utility document from
	the start location to the end location.
	
	If startLocPtr is NULL, the start of the document is used as the start location.
	
	If endLocPtr is NULL, the end of the document is used as the end location.
	
	textHandlePtr is a pointer to your Handle variable.
	
	reservedForFuture must be NULL for now.
	
	flags should be 0 for now.
	
	Example:
		Handle h;
		int result;
		
		if (result = TUFetchText2(TU, NULL, NULL, &h, NULL, 0))	// Fetch all text in document
			return result;
		
		<Deal with the text>
		
		DisposeHandle(h);
	
	Note that the text in the handle h is NOT null terminated and therefore
	is not a C string. You can turn it into a C string as follows:
		length = GetHandleSize(h);
		SetHandleSize(h, length+1);
		if (result = MemError()) {
			DisposeHandle(h);
			return result
		}
		*h[length] = 0;

	Returns 0 if OK, an Igor error code if an error occurs fetching the text
	or the version of Igor that you are running with does not support this callback.
	
	Added in Igor Pro 6.20 but works with any version.
	
	Thread Safety: TUFetchText2 is not thread-safe.
*/
int
TUFetchText2(TUStuffHandle TU, TULocPtr startLocPtr, TULocPtr endLocPtr, Handle* textHandlePtr, void* reservedForFuture, int flags)
{
	if (!CheckRunningInMainThread("TUFetchText2"))
		return NOT_IN_THREADSAFE;
		
	if (igorVersion < 620) {				// Emulate for old versions of Igor
		TULoc saveStartLoc, saveEndLoc;
		TULoc startLoc, endLoc;
		int err;

		TUGetSelLocs(TU, &saveStartLoc, &saveEndLoc);
		
		if (startLocPtr == NULL) {
			startLoc.paragraph = 0;
			startLoc.pos = 0;
		}
		else {
			startLoc = *startLocPtr;
		}
		
		if (endLocPtr == NULL) {
			TUDocInfo di;
			int length;
			
			di.version = TUDOCINFO_VERSION;
			if (err = TUGetDocInfo(TU, &di))
				return err;
			endLoc.paragraph = di.paragraphs-1;
				
			if (err = TUFetchParagraphText(TU, endLoc.paragraph, NULL, &length))
				return err;
			endLoc.pos = length;
		}
		else {
			endLoc = *endLocPtr;
		}

		if (err = TUSetSelLocs(TU, &startLoc, &endLoc, 0))
			return err;
		
		err = TUFetchSelectedText(TU, textHandlePtr, NULL, 0);

		TUSetSelLocs(TU, &saveStartLoc, &saveEndLoc, 0);

		return err;
	}

	return (int)CallBack6(TUFETCHTEXT2, TU, startLocPtr, endLocPtr, textHandlePtr,  reservedForFuture, (void*)flags);
}

/*	TUSetStatusArea(TU, message, eraseFlags, statusAreaWidth)

	If message is not NULL, sets the status message in the text utility document.
	message is a C string. Only the first 127 characters are displayed.
	
	If message is not NULL then eraseFlags determines when the status
	message will be erased. See the TU_ERASE_STATUS #defines in IgorXOP.h.
	
	If statusAreaWidth is >= 0, sets the width of the status area.
	This is in pixels.

	Returns 0 if OK, an Igor error code if the version of Igor that you are running
	with does not support this callback.
	
	Thread Safety: TUSetStatusArea is not thread-safe.
*/
int
TUSetStatusArea(TUStuffHandle TU, const char* message, int eraseFlags, int statusAreaWidth)
{
	if (!CheckRunningInMainThread("TUSetStatusArea"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(TUSETSTATUSAREA, TU, (void*)message,  (void*)eraseFlags, (void*)statusAreaWidth);
}

/*	TUMoveToPreferredPosition(TUStuffHandle TU)

	Moves the window to the preferred position, as determined by the user's notebook
	preferences. Normally, you will call this in response to the MOVE_TO_PREFERRED_POSITION
	message from IGOR.
	
	During the TUMoveToPreferredPosition call, your XOP will receive a GROW message from IGOR.
	On Windows, your window procedure may also receive several messages from the operating system.
	
	Thread Safety: TUMoveToPreferredPosition is not thread-safe.
*/
void
TUMoveToPreferredPosition(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TUMoveToPreferredPosition"))
		return;

	CallBack1(TUMOVE_TO_PREFERRED_POSITION, TU);
}

/*	TUMoveToFullSizePosition(TUStuffHandle TU)

	Moves the window to show all of its content or to fill the screen (Macintosh) or
	MDI frame window (Windows). Normally, you will call this in response to the
	MOVE_TO_FULLSIZE_POSITION message from IGOR.
	
	During the TUMoveToFullSizePosition call, your XOP will receive a GROW message from IGOR.
	On Windows, your window procedure may also receive several messages from the operating system.
	
	Thread Safety: TUMoveToFullSizePosition is not thread-safe.
*/
void
TUMoveToFullSizePosition(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TUMoveToFullSizePosition"))
		return;

	CallBack1(TUMOVE_TO_FULL_POSITION, TU);
}

/*	TURetrieveWindow(TUStuffHandle TU)

	Moves the window, if necessary, to fit entirely within the screen (Macintosh) or
	MDI frame window (Windows). Normally, you will call this in response to the
	RETRIEVE message from IGOR.
	
	During the TURetrieveWindow call, your XOP will receive a GROW message from IGOR.
	On Windows, your window procedure may also receive several messages from the operating system.
	
	Thread Safety: TURetrieveWindow is not thread-safe.
*/
void
TURetrieveWindow(TUStuffHandle TU)
{
	if (!CheckRunningInMainThread("TURetrieveWindow"))
		return;

	CallBack1(TURETRIEVE, TU);
}

/*	HistoryDisplaySelection()

	Scrolls the current selection in the history area into view.
	
	Thread Safety: HistoryDisplaySelection is not thread-safe.
*/
void
HistoryDisplaySelection(void)
{
	if (!CheckRunningInMainThread("HistoryDisplaySelection"))
		return;

	CallBack0(HISTORY_DISPLAYSELECTION);
}

/*	HistoryInsert()

	Inserts the specified text into the history area, replacing the current selection, if any.
	
	If you just want to append text to the history, call XOPNotice or XOPNotice2 instead of
	HistoryInsert.
	
	Except in very rare cases you should not modify the history, except to append to it.
	
	Thread Safety: HistoryInsert is not thread-safe.
*/
void
HistoryInsert(const char* dataPtr, int dataLen)
{
	if (!CheckRunningInMainThread("HistoryInsert"))
		return;

	CallBack2(HISTORY_INSERT, (void*)dataPtr, (void*)dataLen);
}

/*	HistoryDelete()

	Deletes the currently selected text in the history area.
	
	Except in very rare cases you should not modify the history, except to append to it.
	
	Thread Safety: HistoryDelete is not thread-safe.
*/
void
HistoryDelete(void)
{
	if (!CheckRunningInMainThread("HistoryDelete"))
		return;

	CallBack0(HISTORY_DELETE);
}

/*	HistoryLines()

	Returns the number of lines of text in the history area.
	
	Thread Safety: HistoryLines is not thread-safe.
*/
int
HistoryLines(void)
{
	if (!CheckRunningInMainThread("HistoryLines"))
		return 0;

	if (igorVersion < 600)
		return 0;
	
	return (int)CallBack0(HISTORY_LINES);
}

/*	HistoryGetSelLocs()

	Sets *startLocPtr and *endLocPtr to describe the selected text in the history area.

	Returns 0 if OK or an Igor error code.
	
	Thread Safety: HistoryGetSelLocs is not thread-safe.
*/
int
HistoryGetSelLocs(TULocPtr startLocPtr, TULocPtr endLocPtr)
{
	if (!CheckRunningInMainThread("HistoryGetSelLocs"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack2(HISTORY_GETSELLOCS, startLocPtr, endLocPtr);
}

/*	HistorySetSelLocs()

	Sets the selection in the history area.

	If startLocPtr is NULL, the start location is taken to be the start of history area.

	If endLocPtr is NULL, the end location is taken to be the end of history area.
	
	If flags is 1, displays the selection if it is out of view.
	Other bits in flags must be set to 0 as they may be used for other purposes in the future.

	Returns 0 if OK or an Igor error code.
	
	Returns an error if the start or end locations are out of bounds or if the start location is
	after the end location.
	
	Thread Safety: HistorySetSelLocs is not thread-safe.
*/
int
HistorySetSelLocs(TULocPtr startLocPtr, TULocPtr endLocPtr, int flags)
{
	if (!CheckRunningInMainThread("HistorySetSelLocs"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(HISTORY_SETSELLOCS, startLocPtr, endLocPtr, (void*)flags);
}

/*	HistoryFetchParagraphText()

	Like TUFetchParagraphText but it acts on the history area. See TUFetchParagraphText documentation for details.
	
	To get the length of the paragraph text without actually getting the text, pass NULL for textPtrPtr.

	I textPtrPtr is not NULL and *textPtrPtr is not NULL then you must dispose *textPtrPtr using
	DisposePtr when you no longer need it.
	
	Thread Safety: HistoryFetchParagraphText is not thread-safe.
*/
int
HistoryFetchParagraphText(int paragraph,  Ptr* textPtrPtr, int* lengthPtr)
{
	if (textPtrPtr != NULL)
		*textPtrPtr = NULL;

	if (!CheckRunningInMainThread("HistoryFetchParagraphText"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(HISTORY_FETCHPARAGRAPHTEXT, (void*)paragraph, textPtrPtr, lengthPtr);
}

/*	HistoryFetchText(startLocPtr, endLocPtr, textHPtr)

	Returns the history area text from the specified start location to the specified end location.

	If startLocPtr is NULL, the start location is taken to be the start of history area.

	If endLocPtr is NULL, the end location is taken to be the end of history area.
	
	On return, if there is an error, *textHPtr will be NULL. If there is no error, *textHPtr
	will point to a handle containing the text. *textHPtr is not null-terminated. *textHPtr
	belongs to you so dispose it using DisposeHandle when you are finished with it.
	
	Thread Safety: HistoryFetchText is not thread-safe.
*/
int
HistoryFetchText(TULocPtr startLocPtr, TULocPtr endLocPtr, Handle* textHPtr)
{
	*textHPtr = NULL;
	
	if (!CheckRunningInMainThread("HistoryFetchText"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(HISTORY_FETCHTEXT, startLocPtr, endLocPtr, textHPtr);
}
