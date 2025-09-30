/*
 * File:   dsPIC33FJ64GS606TIMER.c
 * Author: Yiou
 *
 * Created on July 7, 2016, 9:55 AM
 */

#include "dsPIC33FJ64GS606TIMER.h"

void initTimer1(int ADCperiod) //working, tested on 7/5/2016
{
    //Scale Timer 1
    //If Fcy = 29.49 MHz, divided by 64, Timer is clocked by 460.78 kHz. Then 1ms is 460;
    T1CONbits.TCKPS = 2;
    PR1 = ADCperiod;
}

void EnableTimer1()
{
    T1CONbits.TON = 1;
}

void initTimer2()
{
    //Scale Timer 2
    //If Fcy = 29.49 MHz, divided by 64, Timer is clocked by 460k kHz. Then 1ms is 460;
    T2CONbits.TCKPS = 2;
    PR2 = IC1timer2Period; //~3 ms
}

void EnableTimer2()
{
    T2CONbits.TON = 1;
}

void initIC1()
{
    IC1CONbits.ICTMR = 1; //select Timer 2 as the source
    IC1CONbits.ICI = 1;   //Interrupt generate every second event (second rising edge)
    IC1CONbits.ICM = 1;   //Mode is to detect every rising and falling edge.
    
    // Enable Capture Interrupt And Timer2
    IPC0bits.IC1IP = 1; // Setup IC1 interrupt priority level
    IFS0bits.IC1IF = 0; // Clear IC1 Interrupt Status Flag
    IEC0bits.IC1IE = 1; // Enable IC1 interrupt
}

void initTimer45(int EndTimeMSB,int EndTimeLSB) //working, tested on 7/5/2016
{
    //If Fcy = 29.49 MHz, divided by 256, Timer is clocked by 115.195 kHz. Then 1ms is 115.195, 230.3906;
    T4CONbits.T32 = 1;
    T4CONbits.TCKPS = 3;
    T5CONbits.TCKPS = 3;
    TMR4 = 0;
    TMR5 = 0;
    PR4 = EndTimeLSB;
    PR5 = EndTimeMSB;
    
    //Set up interrupt
    IFS1bits.T4IF = 0;
    IFS1bits.T5IF = 0;
    IEC1bits.T4IE = 0;  //should not set to 1
    IEC1bits.T5IE = 1;
}

void EnableTimer45()
{
    T4CONbits.TON = 1;
    T5CONbits.TON = 1;
}

void DisableTimer45()
{
    T4CONbits.TON = 0;
    T5CONbits.TON = 0;
    TMR4 = 0;
    TMR5 = 0;
}