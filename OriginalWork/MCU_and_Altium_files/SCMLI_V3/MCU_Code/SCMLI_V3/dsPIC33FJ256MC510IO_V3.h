/* 
 * File:   dsPIC33FJ64GS606IO.h
 * Author: Yiou He
 * Comments: header file of dsPIC33FJ64GS606 IO file
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
// more than once.  
#ifndef XC_HEADER_IO_H
#define	XC_HEADER_IO_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ256MC510.h"
#include "stdio.h"
#include "string.h"

void initIO();
void SetModules(int *state_indices, int index);
void TestOut(int sel);
int TestRead(int sel_int);
void DirectSet(int sel, int val);
void QuickSet(int index);
#endif	/* XC_HEADER_IO_H */