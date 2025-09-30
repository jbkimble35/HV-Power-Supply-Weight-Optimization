/* 
 * File:   dsPIC33FJ64GS606SPI.h
 * Author: Yiou He
 * Comments: header file of dsPIC33FJ64GS606 SPI unit
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
// more than once.  
#ifndef XC_HEADER_SPI_H
#define	XC_HEADER_SPI_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ64GS606.h"
#include "stdio.h"
#include "string.h"
#include "main_PWMoutputDefinitions.h"

#define WR_ENABLE 6
#define WR_DISABLE 4
#define RD_STATUS 5
#define WR_STATUS 1
#define READ 3
#define WRITE 2
#define RDID 0x9F
#define RAMCS LATGbits.LATG9
#define RAM_enable() RAMCS = 0
#define RAM_disable() RAMCS = 1;

void initSPI2();
void EnableSPI2();
void DisableSPI2();
void writeSPI2 (unsigned int in);
unsigned int readSPI2();
void write_enable(void);
void write_disable(void);
void writeByte(unsigned int in,int long loc);
unsigned int readByte(int long loc);
int long ReadID();
int long SendThroughSPI(int Recorder,double InputVoltage,double InputCurrent,double OutputVoltage,double Temperature,int long loc);

#endif	/* XC_HEADER_SPI_H */