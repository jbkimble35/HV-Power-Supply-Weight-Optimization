/*
 * File:   dsPIC33FJ64GS606IO.c
 * Author: Yiou
 *
 * Created on July 7, 2016, 9:53 AM
 */

#include "dsPIC33FJ64GS606IO.h"

void initIO()  //working, tested on 7/4/2016
{
    //Initialize the IO port
    //Use PORT to read, use LAT to write
    //PortB (Pin 0-4 (ANALOG), 6-7 (Enables), 8 (BT_DSR), 9-13, 15 as output, 14 (BT_DTR) as input)
    LATB = 0x0000;    //Set the latch to zeros for all
    TRISB = 0x401F; //0b0100000000011111;
    //PORTC all output (RC12 & 14 as PGED and C)
    LATC = 0x0000;
    TRISC = 0x0000;
    //PORTD all output (NOT USED)
    LATD = 0x0000;  //RD4 - 7: 200VON, TBD, HV, Logic
    TRISD = 0x0000;
    //PORTE all output (RE0 - 3 as PWM output)
    LATE = 0x0000;
    TRISE = 0x0000;
    //PORTF (RF4 - 5 are RX and TX, RF2 - 3 are RX and TX, RF6 is RCStop input)
    TRISF = 0x0014; //0b0000000000010100
    LATF = 0x0000;
    //PORTG (RG6-9 are SCK, SI, SO, CS, SI (RG7) is input)
    LATG = 0x0000;
    TRISG = 0x0080;  //0b0000000010000000
    //Analog and comparator input setting
    //All analog and comparator pins are default as analog input, need to be set to 1 when it's used as digital pins
    //Using AN0 - 4 as analog and comparator input, using AN 6-8, 10, 11, 14 as digital IO, not using AN5, 9, 15
    ADPCFG = 0xFFE0;
    //Disable comparator pairs
}

void DischargeCaps() //Start to discharge the input caps
{
    //LOGIC IS ON
    LATBbits.LATB6 = 0;
    LATBbits.LATB7 = 0;
    LATDbits.LATD4 = 0;
    LATDbits.LATD5 = 0;
    LATDbits.LATD6 = 0;
    LATDbits.LATD7 = 1;
}

void Precharge()      //Precharge the input caps
{
    //200V is On
    LATBbits.LATB6 = 1;
    LATBbits.LATB7 = 0;
    LATDbits.LATD4 = 1;
    LATDbits.LATD5 = 0;
    LATDbits.LATD6 = 0;
    LATDbits.LATD7 = 0;
}

void NormalOperation()  //Connect 200V to the switches
{
    //200V is on
    LATBbits.LATB6 = 0;
    LATBbits.LATB7 = 1;
    LATDbits.LATD4 = 1;
    LATDbits.LATD5 = 0;
    LATDbits.LATD6 = 0;
    LATDbits.LATD7 = 0;
}

void HVon()
{
    //HV is on
    LATDbits.LATD4 = 0;
    LATDbits.LATD5 = 0;
    LATDbits.LATD6 = 1;
    LATDbits.LATD7 = 0;
}