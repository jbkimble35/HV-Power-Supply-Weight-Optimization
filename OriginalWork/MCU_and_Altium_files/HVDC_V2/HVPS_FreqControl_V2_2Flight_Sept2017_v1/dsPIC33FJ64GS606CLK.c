/*
 * File:   dsPIC33FJ64GS606CLK.c
 * Author: Yiou
 *
 * Created on July 7, 2016, 9:55 AM
 */

#include "dsPIC33FJ64GS606CLK.h"

void initCLK() //working, tested on 7/5/2016
{
    //Tune the OSC
    OSCTUNbits.TUN = 0b111011; //Tune the frequency to lower
    /* Configure Oscillator to operate the device at 29.84 Mhz
	   Fosc= Fin*M/(N1*N2), Fcy=Fosc/2
 	   Fosc= 7.37*(64)/(4*2) = Mhz for Fosc, Fcy = 29.84 Mhz
       Under this condition, Icc = 61 mA) */
	/* Configure PLL prescaler, PLL postscaler, PLL divisor */
	
    PLLFBD = 62; 		    /* M = PLLFBD + 2, working */
	CLKDIVbits.PLLPOST = 1;   /* N1 = 4, working, if PLLPOST = 0, N1 = 2*/
	CLKDIVbits.PLLPRE = 0;    /* N2 = 2, working */
    //Check FRCDIV, working	CLKDIVbits.FRCDIV = 1;  /* Divide FRC with 2*/
    while(OSCCONbits.LOCK != 1);			/* Wait for Pll to Lock */
    //OSCCONbits.NOSC = 
    //Reference CLK output (OSC2 should be set as general IO to enable this function, working)
    //REFOCONbits.ROSEL = 0;  // Reference Clock Source Select System Clock
    //REFOCONbits.RODIV = 15;  //Divide the main clock by 32768, Fosc is 80MHz, Reference CLK should be 2.41 kHz 
    //REFOCONbits.ROON = 1;   //Reference Oscillator is enabled at pin 40	
    
    /* Now setup the ADC and PWM clock for 120MHz
	   ((FRC * 16) / APSTSCLR ) = (7.37 * 16) / 1 = ~ 120MHz*/
	ACLKCONbits.FRCSEL = 1;				/* FRC provides input for Auxiliary PLL (x16) */
	ACLKCONbits.SELACLK = 1;		    /* Auxiliary Oscillator provides clock source for PWM & ADC */
	ACLKCONbits.APSTSCLR = 7;			/* Divide Auxiliary clock by 1 */
	ACLKCONbits.ENAPLL = 1;				/* Enable Auxiliary PLL */
    while(ACLKCONbits.APLLCK != 1);		/* Wait for Auxiliary Pll to Lock */
}