/*
 * File:   dsPIC33FJ64GS606SPI.c
 * Author: Yiou
 *
 * Created on July 7, 2016, 9:44 AM
 */

#include "dsPIC33FJ64GS606SPI.h"

//init SPI 2 Module for external nonvolatile RAM 
void initSPI2()
{
    IFS2bits.SPI2IF = 0; // Clear the Interrupt flag
    IEC2bits.SPI2IE = 0; // Disable the interrupt
    
    SPI2STAT = 0; //Clear status register
    SPI2CON1 = 0; //Clear control register 1
    SPI2CON2 = 0; //Clear control register 2
    
    // SPI2CON1 Register Settings
    SPI2CON1bits.DISSCK = 0; // Internal serial clock is enabled
    SPI2CON1bits.DISSDO = 0; // SDOx pin is controlled by the module
    SPI2CON1bits.MODE16 = 0; // Communication is byte-wide (8 bits)
    SPI2CON1bits.MSTEN = 1; // Master mode enabled
    SPI2CON2bits.FRMEN = 0; //Frame mode
    SPI2CON2bits.FRMPOL = 0;
    SPI2CON1bits.SMP = 0; // Input data is sampled at the middle of data output time
    SPI2CON1bits.CKE = 1; // Serial output data changes on transition from
    SPI2CON1bits.PPRE = 3; //Postscale clock to 1:1
    SPI2CON1bits.SPRE = 5;  //Prescale clock to 1:1
    // Idle clock state to active clock state
    SPI2CON1bits.CKP = 0; // Idle state for clock is a low level; 

    // Interrupt Controller Settings
    //C1FEN1bits.FLTEN11 = 0;
    IFS2bits.SPI2IF = 0;  // Clear the Interrupt flag
    IEC2bits.SPI2IE = 0;  // Disable the interrupt
}

void EnableSPI2()
{
    SPI2STATbits.SPIEN = 1;
}

void DisableSPI2()
{
    SPI2STATbits.SPIEN = 0;
}

//send one byte of data and receive one back at the same time
 void writeSPI2(unsigned int in)
 {
    int i = 0;
    SPI2STATbits.SPIROV = 0;
    SPI2BUF = in;					// write to buffer for TX
    while(SPI2STATbits.SPITBF);	// wait for transfer to complete
    //in = SPI2BUF;    				// avoid overflow when reading
    for (i = 0; i<= 1;i++);
 } // writeSPI

 void write_enable(void)
 {
     RAM_enable();
     writeSPI2(WR_ENABLE);
     RAM_disable();
 }
 
 void write_disable(void)
 {
    RAM_enable();
    writeSPI2(WR_DISABLE);
    RAM_disable();
 }

 void writeByte(unsigned int in,int long loc)
 {
     RAM_enable();
     writeSPI2(WR_ENABLE);
     RAM_disable();  //     write_enable();
     RAM_enable();
     writeSPI2(WRITE);
     writeSPI2(loc >> 16);
     writeSPI2(loc >> 8);
     writeSPI2(loc);
     writeSPI2(in);
     RAM_disable();
 }
 
 unsigned int readByte(int long loc)
 { 
     int i = 0;
     unsigned int r;
     RAM_enable();
     writeSPI2(READ);
     writeSPI2(loc >> 16);
     writeSPI2(loc >> 8);
     writeSPI2(loc);
     r = readSPI2();
     for (i = 0; i<= 1;i++);
     RAM_disable();
     return r;
 }

   unsigned int readSPI2()
{
    //int i = 0;
    unsigned int out = -1;
    SPI2STATbits.SPIROV = 0;
    SPI2BUF = 0x13;                  // initiate bus cycle 
    while(!SPI2STATbits.SPIRBF);
    /*Check for Receive buffer full status bit of status register*/
    if (SPI2STATbits.SPIRBF)
    {
        //SPI2STATbits.SPIROV = 0;
        //for (i = 0; i<= 0;i++);
        out = SPI2BUF;    /* return byte read */
    }

    return out;                  		/* RBF bit is not set return error*/
}
  
int long ReadID()
{
    int long r = 0;
    RAM_enable();
    writeSPI2(RDID);
    r = readSPI2();
    RAM_disable();
    return r;
}

int long SendThroughSPI(int Recorder, double InputVoltage,double InputCurrent,double OutputVoltage,double Temperature,int long loc)
{
    writeByte(Recorder,loc);
    loc++;
    writeByte(((unsigned int)InputVoltage) >> 2,loc);  //200 V corresponds to 626
    loc++;
    writeByte(((unsigned int)InputCurrent) >> 2,loc);
    loc++;
    writeByte(((unsigned int)OutputVoltage) >> 2,loc); //40 kV corresponds to 682
    loc++;
    writeByte(((unsigned int)Temperature) >> 1,loc); //usually around 300, maximum with hot water is 375
    loc++;
    if (loc >= MaxLoc) loc = MaxLoc;
    return loc;
}