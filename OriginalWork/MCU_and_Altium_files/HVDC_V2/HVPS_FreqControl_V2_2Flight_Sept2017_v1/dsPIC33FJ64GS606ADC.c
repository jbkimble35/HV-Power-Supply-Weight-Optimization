/*
 * File:   dsPIC33FJ64GS606ADC.c
 * Author: Yiou
 *
 * Created on July 7, 2016, 9:41 AM
 */

#include "dsPIC33FJ64GS606ADC.h"

void initADC(int period)  
{
    //We will use AN 0 -6, pair 0 - 3
    ADCONbits.SLOWCLK = 1;   //ADC is clocked by auxiliary PLL.
    ADCONbits.FORM = 0;	     //Integer data format
    ADCONbits.EIE = 0;		 //Interrupt generated after the first conversion is completed
    ADCONbits.ORDER = 0;	 //Convert even channel first
    ADCONbits.SEQSAMP = 0;	 //Shared S&H circuit and dedicated simultaneous sampling
    ADCONbits.ADCS = 5;		 //ADC clock = FADC/8 = 120MHz / 8 = 15MHz

    //Set PWM as the trigger source of ADC (changed to Timer 1)
   	ADCPC0bits.TRGSRC0 = 12; 		//ADC Pair 0 triggered by Timer 1
    ADCPC0bits.TRGSRC1 = 12; 		//ADC Pair 1 triggered by Timer 1
    ADCPC1bits.TRGSRC2 = 12; 		//ADC Pair 2 triggered by Timer 1
    //ADCPC1bits.TRGSRC3 = 12; 		//ADC Pair 3 triggered by Timer 1
    
    //Clear pair ready bit
    ADSTATbits.P0RDY = 0; 			//Clear all pairs data ready bit
    ADSTATbits.P1RDY = 0; 			//Clear all pairs data ready bit
    ADSTATbits.P2RDY = 0; 			//Clear all pairs data ready bit
    //ADSTATbits.P3RDY = 0; 			//Clear all pairs data ready bit

    //Set interrupt
    //ADCPC0bits.IRQEN0 = 1;		    //Enable ADC Interrupt pair 0 
    IFS6bits.ADCP0IF = 0;			//Clear ADC Pair 0 interrupt flag
   	IEC6bits.ADCP0IE = 1;			//Enable the ADC Pair 0 interrupt
    IFS6bits.ADCP1IF = 0;			//Clear ADC Pair 1 interrupt flag
   	IEC6bits.ADCP1IE = 1;			//Enable the ADC Pair 1 interrupt
    IFS7bits.ADCP2IF = 0;			//Clear ADC Pair 2 interrupt flag
   	IEC7bits.ADCP2IE = 1;			//Enable the ADC Pair 2 interrupt
    //IFS7bits.ADCP3IF = 0;			//Clear ADC Pair 3 interrupt flag
   	//IEC7bits.ADCP3IE = 1;			//Enable the ADC Pair 3 interrupt

    //Set the PWM trigger
    //TRGCON1bits.TRGDIV = 15;  //Set the trigger event to divide by 16
 	//TRGCON1bits.TRGSTRT = 0;  // enable Trigger generated after 0 PWM cycles
    //TRIG1 = period; //Trigger every time when it reaches to the period.
}

void EnableADC()
{
    ADCONbits.ADON = 1;
}

void DisableADC()
{
    ADCONbits.ADON = 0;
}