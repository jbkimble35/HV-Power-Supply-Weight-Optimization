/* 
 * File: dsPIC33FJ64GS606UART.h
 * Author: Yiou He
 * Comments: header file for the UART function
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
#ifndef XC_HEADER_UART_H
#define	XC_HEADER_UART_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ64GS606.h"
#include "stdio.h"
#include "string.h"
#include "main_PWMoutputDefinitions.h"

void initUART1();
void initUART2();
void EnableUART1();  //StartUART1
void EnableUART2();   //StartUART2
void DisableUART1();
void WriteScreen(char s[50]); //Write texts on the screen
void WriteBluetooth(char s[50]); //Write texts to the bluetooth
void WriteScreenADCdata(int,double,double,double,double,double,double,double,double);

#endif	/* XC_HEADER_UART_H */