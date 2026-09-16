
// starmap2.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'pch.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// Cstarmap2App:
// See starmap2.cpp for the implementation of this class

class Cstarmap2App : public CWinApp
{
public:
	Cstarmap2App();

// Overrides
public:
	virtual BOOL InitInstance();

// Implementation
	DECLARE_MESSAGE_MAP()
};
