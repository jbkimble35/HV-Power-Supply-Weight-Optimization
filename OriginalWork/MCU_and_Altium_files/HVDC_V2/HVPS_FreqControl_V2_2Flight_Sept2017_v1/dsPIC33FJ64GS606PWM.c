/*
 * File:   dsPIC33FJ64GS606PWM.c
 * Author: Yiou
 *
 * Created on July 7, 2016, 9:48 AM
 */


#include "dsPIC33FJ64GS606PWM.h"

void initPWM(int period, int duty, int deadtime)
{
    PTCON = 0;     //Disable PWM module
    PTCON2 = 0;    //Set the PWM resolution to max 1.04 ns
    PTPER = period;  //Set the primary time base (PWM period) to 500 kHz, 1/500000/1.04ns)
    MDC = duty;     //Set the PWM duty cycle to 50%
    //Initiate PWM1
    PWMCON1 = 0x0000; //Individual PDC control duty cycle, positive dead time
    PDC1 = duty;       //Set the secondary PWM duty cycle
    PHASE1 = 0;
    DTR1 = deadtime;        //PWMH deadtime is DTR1*1.04 ns
    ALTDTR1 = deadtime;     //PWML deadtime is ALTDTR1*1.04 ns
    IOCON1 = 0x0341;       //Override PWM, not controlled by PWM module
    //Initiate PWM2
    PWMCON2 = 0x0000;     //Individual PDC control duty cycle
    PDC2 = period;
    PHASE2 = 0;
    DTR2 = deadtime;      //PWMH deadtime is DTR2*1.04 ns
    ALTDTR2 = deadtime;   //PWML deadtime is ALTDTR2*1.04 ns
    IOCON2 = 0x0341;      //Override PWM, not controlled by PWM module
    /*
     ;GPIO module controls PWMxH pin
	 ;GPIO module controls PWMxL pin
	 ;complementary PWM output mode
	 ;overide enable for PWMxH
	 ;overide enable for PWMxL
     ;override date 01
	 ;overide syncronized to PWM time base
     */
}
void EnablePWM()
{
    /*
    TRISEbits.TRISE1 = 1; // Configure PWM1H/RE1 as digital input
    TRISEbits.TRISE0 = 1; // Configure PWM1L/RE0 as digital input
    TRISEbits.TRISE3 = 1; // Configure PWM2H/RE3 as digital input
    TRISEbits.TRISE2 = 1; // Configure PWM2L/RE2 as digital input
    // Ensure output is in safe state using pull-up or pull-down resistors
    // Did not add pull-down resistors in the version, seems fine?
    
    IOCON1bits.PENH = 0; // Assign pin ownership of PWM1H/RE1 to GPIO module
    IOCON1bits.PENL = 0; // Assign pin ownership of PWM1L/RE0 to GPIO module
    IOCON1bits.OVRDAT = 1; // Configure override state of the PWM outputs to
    IOCON2bits.PENH = 0; // Assign pin ownership of PWM1H/RE1 to GPIO module
    IOCON2bits.PENL = 0; // Assign pin ownership of PWM1L/RE0 to GPIO module
    IOCON2bits.OVRDAT = 1; // Configure override state of the PWM outputs to
    // desired safe state.
    IOCON1bits.OVRENH = 1; // Override PWM1H output
    IOCON1bits.OVRENL = 1; // Override PWM1L output
    IOCON2bits.OVRENH = 1; // Override PWM1H output
    IOCON2bits.OVRENL = 1; // Override PWM1L output
    PTCONbits.PTEN = 1; // Enable PWM module
    */

    //Remove override
    IOCON1bits.OVRENH = 0; // Remove override for PWM1H output
    IOCON1bits.OVRENL = 0; // Remove override for PWM1L output
    IOCON2bits.OVRENH = 0; // Remove override for PWM1H output
    IOCON2bits.OVRENL = 0; // Remove override for PWM1L output
     
    //delay(100); //for(int i = 0;i<= 10; i++); // Introduce a delay greater than one full PWM cycle
    IOCON1bits.PENH = 1; // Assign pin ownership of PWM1H/RE1 to PWM module
    IOCON1bits.PENL = 1; // Assign pin ownership of PWM1L/RE0 to PWM module
    IOCON2bits.PENH = 1; // Assign pin ownership of PWM1H/RE3 to PWM module
    IOCON2bits.PENL = 1; // Assign pin ownership of PWM1L/RE2 to PWM module
    
}

void DisablePWM()
{
    //Disable PWM
    PTCONbits.PTEN = 0;
    // Have to add the similar procedure to make sure the switches are turned off to the correct states.
    TRISEbits.TRISE1 = 0; // Configure PWM1H/RE1 as digital output
    TRISEbits.TRISE0 = 0; // Configure PWM1L/RE0 as digital output
    TRISEbits.TRISE3 = 0; // Configure PWM2H/RE3 as digital output
    TRISEbits.TRISE2 = 0; // Configure PWM2L/RE2 as digital output
    
    LATEbits.LATE1 = 0; //PWM1H = 0
    LATEbits.LATE0 = 1; //PWM1L = 1
    LATEbits.LATE3 = 1; //PWM2H = 0
    LATEbits.LATE2 = 0; //PWM2L = 1
    /*
    TRISEbits.TRISE1 = 0; // Configure PWM1H/RE1 as digital input
    TRISEbits.TRISE0 = 0; // Configure PWM1L/RE0 as digital input
    TRISEbits.TRISE3 = 0; // Configure PWM2H/RE3 as digital input
    TRISEbits.TRISE2 = 0; // Configure PWM2L/RE2 as digital input
    // Ensure output is in safe state using pull-up or pull-down resistors
    // Did not add pull-down resistors in the version, seems fine?
    
    IOCON1bits.PENH = 0; // Assign pin ownership of PWM1H/RE1 to GPIO module
    IOCON1bits.PENL = 0; // Assign pin ownership of PWM1L/RE0 to GPIO module
    IOCON1bits.OVRDAT = 1; // Configure override state of the PWM outputs to
    IOCON2bits.PENH = 0; // Assign pin ownership of PWM2H/RE1 to GPIO module
    IOCON2bits.PENL = 0; // Assign pin ownership of PWM2L/RE0 to GPIO module
    IOCON2bits.OVRDAT = 2; // Configure override state of the PWM outputs to
    
    PTCONbits.PTEN = 0; // Disable PWM module
    */
}

/*void PrechargePWM()
{
    TRISEbits.TRISE1 = 0; // Configure PWM1H/RE1 as digital input
    TRISEbits.TRISE0 = 0; // Configure PWM1L/RE0 as digital input
    TRISEbits.TRISE3 = 0; // Configure PWM2H/RE3 as digital input
    TRISEbits.TRISE2 = 0; // Configure PWM2L/RE2 as digital input
    
    //Control PWM output with GPIO
    LATEbits.LATE1 = 0; //PWM1H = 0
    LATEbits.LATE0 = 1; //PWM1L = 1
    LATEbits.LATE3 = 1; //PWM2H = 0
    LATEbits.LATE2 = 0; //PWM2L = 1
}*/

void PreparePWM()
{
    TRISEbits.TRISE1 = 0; // Configure PWM1H/RE1 as digital input
    TRISEbits.TRISE0 = 0; // Configure PWM1L/RE0 as digital input
    TRISEbits.TRISE3 = 0; // Configure PWM2H/RE3 as digital input
    TRISEbits.TRISE2 = 0; // Configure PWM2L/RE2 as digital input
    // Ensure output is in safe state using pull-up or pull-down resistors
    // Did not add pull-down resistors in the version, seems fine?
    
    IOCON1bits.PENH = 0; // Assign pin ownership of PWM1H/RE1 to GPIO module
    IOCON1bits.PENL = 0; // Assign pin ownership of PWM1L/RE0 to GPIO module
    IOCON1bits.OVRDAT = 1; // Configure override state of the PWM outputs to
    IOCON2bits.PENH = 0; // Assign pin ownership of PWM2H/RE1 to GPIO module
    IOCON2bits.PENL = 0; // Assign pin ownership of PWM2L/RE0 to GPIO module
    IOCON2bits.OVRDAT = 2; // Configure override state of the PWM outputs to
    // desired safe state.
    IOCON1bits.OVRENH = 1; // Override PWM1H output
    IOCON1bits.OVRENL = 1; // Override PWM1L output
    IOCON2bits.OVRENH = 1; // Override PWM1H output
    IOCON2bits.OVRENL = 1; // Override PWM1L output
    
    PTCONbits.PTEN = 1; // Enable PWM module
}

void ChangePD(int period, int duty, int deadtime)
{
    PTPER = period;  //Set the primary time base
    MDC = duty;     //Set the PWM duty cycle to 50%
    PDC1 = duty;       //Set the secondary PWM duty cycle
    PDC2 = duty;       //Set the secondary PWM duty cycle
    
    DTR1 = deadtime;        //PWMH deadtime is DTR1*1.04 ns
    ALTDTR1 = deadtime;     //PWML deadtime is ALTDTR1*1.04 ns
    DTR2 = deadtime;      //PWMH deadtime is DTR2*1.04 ns
    ALTDTR2 = deadtime;   //PWML deadtime is ALTDTR2*1.04 ns
}
//
void initChangePD(int period, int duty, int deadtime)
{
    PTPER = period;  //Set the primary time base
    MDC = duty;     //Set the PWM duty cycle to 50%
    PDC1 = duty;       //Set the secondary PWM duty cycle
    PDC2 = period;       //Set the secondary PWM duty cycle
    
    DTR1 = deadtime;        //PWMH deadtime is DTR1*1.04 ns
    ALTDTR1 = deadtime;     //PWML deadtime is ALTDTR1*1.04 ns
    DTR2 = deadtime;      //PWMH deadtime is DTR2*1.04 ns
    ALTDTR2 = deadtime;   //PWML deadtime is ALTDTR2*1.04 ns
}