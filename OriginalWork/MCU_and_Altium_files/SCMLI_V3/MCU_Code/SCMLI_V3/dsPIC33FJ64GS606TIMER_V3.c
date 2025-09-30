/*
 * File:   dsPIC33FJ64GS606TIMER_test.c
 * Author: Suzanne
 *
 * Created on November 28, 2018, 2:42 PM
 */

#include "dsPIC33FJ256MC510TIMER_V3.h"

void initTimer1(int SwitchPeriod) //working, tested on 7/5/2016
{
    //Scale Timer 1
    //If Fcy = 29.49 MHz, divided by 64, Timer is clocked by 460.78 kHz. Then .01ms is 5;
    T1CONbits.TCKPS = 2;
    TMR1 = 0;
    PR1 = SwitchPeriod;
    
    IFS0bits.T1IF = 0;
    IEC0bits.T1IE = 1;
    
}

void EnableTimer1()
{
    T1CONbits.TON = 1;
}

void DisableTimer1()
{
    T1CONbits.TON = 0;
}

void initTimer2(int StartPeriod) //working, tested on 7/5/2016
{
    //Scale Timer 1
    //If Fcy = 29.49 MHz, divided by 256, Timer is clocked by 11.519 kHz. Then 11.5 is 1ms;
    T2CONbits.TCKPS = 3;
    TMR2 = 0;
    PR2 = StartPeriod;
    
    IFS0bits.T2IF = 0;
    IEC0bits.T2IE = 1;
}

void EnableTimer2()
{
    T2CONbits.TON = 1;
}

void DisableTimer2()
{
    T2CONbits.TON = 0;
}
