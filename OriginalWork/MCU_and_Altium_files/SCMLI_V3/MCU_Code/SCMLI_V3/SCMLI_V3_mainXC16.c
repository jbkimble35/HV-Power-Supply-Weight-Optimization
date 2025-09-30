/*
 * File:   test1_6332_mainXC16.c
 * Author: Suzanne
 *
 * Created on November 27, 2018, 9:18 PM
 */


#include "xc.h"

// Config settings
// POSCMOD = HS, FNOSC = PRIPLL, FWDTEN = OFF
// PLLIDIV = DIV_2, PLLMUL = MUL_16
// PBDIV = 8 (default)

#include "xc.h"
#include "p33FJ256GP510.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "dsPIC33FJ256GP510UART_V3.h"
#include "dsPIC33FJ256GP510CLK_V3.h"
#include "dsPIC33FJ256GP510TIMER_V3.h"
#include "dsPIC33FJ256GP510IO_V3.h"
#include "main_SCMLI_V3.h"

/* Configuration Bit Settings */
// <editor-fold defaultstate="collapsed" desc="Configuration Bits Setup">
#pragma config BWRP = WRPROTECT_OFF
/* ----- FBS (0xf80000) --------------------------------------------------
 *   Boot Segment Write Protect:
 *      BWRP_WRPROTECT_ON    Boot Segment is write protected
 *      BWRP_WRPROTECT_OFF   Boot Segment may be written
 * */
#pragma config GWRP = OFF
/* ----- FGS (0xf80004) --------------------------------------------------
 *   General Code Segment Write Protect:
 *      GWRP_ON              General Segment is write protected
 *      GWRP_OFF             General Segment may be written
 * */
#pragma config FNOSC = FRCPLL
#pragma config IESO = OFF
/* ----- FOSCSEL (0xf80006) --------------------------------------------------
;   Oscillator Source Selection:
;     FNOSC_FRC            Internal Fast RC (FRC), 7.37 MHz
;     FNOSC_FRCPLL         Internal Fast RC with PLL (FRCPLL)
;     FNOSC_PRI            Primary Oscillator (XT, HS, EC)
;     FNOSC_PRIPLL         Primary Oscillator (XT, HS, EC) with PLL
;     FNOSC_SOSC           Secondary Oscillator (SOSC)
;     FNOSC_LPRC           Low-Power RC Oscillator (LPRC)
;     FNOSC_FRCDIV16       Internal Fast RC (FRC) Oscillator with divide-by-16
;     FNOSC_FRCDIVN        Internal Fast RC (FRC) Oscillator with postscaler
;
;   Internal External Switch Over Mode:
;     IESO_OFF             Start up with user-selected oscillator source
;     IESO_ON              Start up device with FRC, then switch to user-selected oscillator source
*/
#pragma config POSCMD = NONE
#pragma config OSCIOFNC = OFF
#pragma config FCKSM = CSECME
/* ----- FOSC (0xf80008) --------------------------------------------------
;
;   Primary Oscillator Source:
;     POSCMD_EC            EC (External Clock) Mode
;     POSCMD_XT            XT Crystal Oscillator Mode
;     POSCMD_HS            HS Crystal Oscillator Mode
;     POSCMD_NONE          Primary Oscillator disabled
;
;   OSC2 Pin Function:
;     OSCIOFNC_ON          OSC2 is general purpose digital I/O pin
;     OSCIOFNC_OFF         OSC2 is clock output, it will output Fcy but not function as REFCLKO.
;
;   Clock Switching Mode bits:
;     FCKSM_CSECME         Both Clock switching and Fail-safe Clock Monitor are enabled
;     FCKSM_CSECMD         Clock switching is enabled,Fail-safe Clock Monitor is disabled
;     FCKSM_CSDCMD         Both Clock switching and Fail-safe Clock Monitor are disabled
*/
#pragma config FWDTEN = OFF
/*;   Watchdog Timer Enable:
;     FWDTEN_OFF           Watchdog timer enabled/disabled by user software
;     FWDTEN_ON            Watchdog timer always enabled
*/
#pragma config FPWRT = PWR128
/*;----- FPOR (0xf8000c) --------------------------------------------------
;
;   POR Timer Value:
;     FPWRT_PWR1           Disabled
;     FPWRT_PWR2           2ms
;     FPWRT_PWR4           4ms
;     FPWRT_PWR8           8ms
;     FPWRT_PWR16          16ms
;     FPWRT_PWR32          32ms
;     FPWRT_PWR64          64ms
;     FPWRT_PWR128         128ms
;
;   Enable Alternate SS1 pin bit:
;     ALTSS1_OFF           SS1 is selected as the I/O pin for SPI1
;     ALTSS1_ON            SS1A is selected as the I/O pin for SPI1
;
;   Enable Alternate QEI1 pin bit:
;     ALTQIO_ON            QEA1A, QEB1A, and INDX1A are selected as inputs to QEI1
;     ALTQIO_OFF           QEA1, QEB1, INDX1 are selected as inputs to QEI1
 * */
#pragma config ICS = PGD2
#pragma config JTAGEN = OFF
/*;----- FICD (0xf8000e) --------------------------------------------------
;
;   Comm Channel Select:
;     ICS_NONE             Reserved, do not use
;     ICS_PGD3             Communicate on PGC3/EMUC3 and PGD3/EMUD3
;     ICS_PGD2             Communicate on PGC2/EMUC2 and PGD2/EMUD2
;     ICS_PGD1             Communicate on PGC1/EMUC1 and PGD1/EMUD1
;
;   JTAG Port Enable:
;     JTAGEN_OFF           JTAG is disabled
;     JTAGEN_ON            JTAG is enabled
 */
//#pragma config HYST0 = HYST0
//#pragma config HYST1 = HYST0
/*;----- FCMP (0xf80010) --------------------------------------------------
;
;   Even Comparator Hysteresis Select:
;     HYST0_HYST0          No Hysteresis
;     HYST0_HYST15         15 mV Hysteresis
;     HYST0_HYST30         30 mV Hysteresis
;     HYST0_HYST45         45 mV Hysteresis
;
;   Comparator Hysteresis Polarity (for even numbered comparators):
;     CMPPOL0_POL_RISE     Hysteresis is applied to rising edge
;     CMPPOL0_POL_FALL     Hysteresis is applied to falling edge
;
;   Odd Comparator Hysteresis Select:
;     HYST1_HYST0          No Hysteresis
;     HYST1_HYST15         15 mV Hysteresis
;     HYST1_HYST30         30 mV Hysteresis
;     HYST1_HYST45         45 mV Hysteresis
;
;   Comparator Hysteresis Polarity (for odd numbered comparators):
;     CMPPOL1_POL_RISE     Hysteresis is applied to rising edge
;     CMPPOL1_POL_FALL     Hysteresis is applied to falling edge
*/
// </editor-fold>

#define SHORT_DELAY (50*8)
#define LONG_DELAY (400*8)

void TestStart();
void TestTimer1();
void TestTimer2();
void ChargeFly();
void SimpleConfig(int time, int index);
void DirectConfig(int time, int sel, int val);

const int startup_dt[3] = { 10, 100, 100 }; //1 count is .085ms

const int PeriodCounts = 460; //roughly 0.5ms half-period
const int NumModules = 2;

int FlagTimer1End = 0;
int FlagTimer2End = 0;

int main(void) {
    
    WriteScreen("Starting");
    
    //initialize the IO
    initIO();
    
    //should there be some pre-startup checks here? probably
    //TestTimer1();
    
    int charge_con[11] =  {1,1,1,1, 1,1,1,1, 1,1,1};
    int con_0V[11] =      {2,2,2,2, 2,2,2,2, 2,2,2};
    int charge_boot[11] = {5,5,5,5, 5,5,5,5, 5,5,5};
    int neg[11] =         {4,4,4,4, 4,4,4,4, 4,4,4};
    int all0[11] =        {0,0,0,0, 0,0,0,0, 0,0,0};
    int charge1[11] =     {1,1,1,0, 0,0,0,0, 0,0,0}; //ADE 
    int boot1[11] =       {5,5,5,0, 0,0,0,0, 0,0,0}; //BD
    int charge2[11] =     {1,1,1,1, 0,0,0,0, 0,0,0}; //ADE 
    int boot2[11] =       {5,5,5,5, 0,0,0,0, 0,0,0}; //BD
    int charge3[11] =     {1,1,1,1, 1,0,0,0, 0,0,0}; //ADE 
    int boot3[11] =       {5,5,5,5, 5,0,0,0, 0,0,0}; //BD
    
    //SimpleConfig(5, charge_boot); 
    SetModules(boot2, 0);
    SetModules(charge2, 1); //BD on, high output 
    
    while (1)
    {
        SimpleConfig(1, 0); //all 5s (boot) BD
        SimpleConfig(2, 1); //all 1s (gen1) ADE
    } 
    
    WriteScreen("Done");
    
    return 0;
}

void TestStart() {
    //SetState(0);
}

void TestTimer1() {
    while(1) {
        initTimer1(startup_dt[1]);
        EnableTimer1();
        TestOut(0);
        while (FlagTimer1End == 0); //program is getting hella caught here 
        FlagTimer1End = 0;
        DisableTimer1();
        
        //initTimer1(10);
        initTimer1(startup_dt[2]);
        EnableTimer1();
        TestOut(1);
        while (FlagTimer1End == 0);
        FlagTimer1End = 0;
        DisableTimer1();
    }
}

void TestTimer2() {
    while(1) {
        initTimer2(startup_dt[1]);
        EnableTimer2();
        TestOut(0);
        while (FlagTimer2End == 0);
        FlagTimer2End = 0;
        DisableTimer2();
        
        initTimer2(startup_dt[2]);
        EnableTimer2();
        TestOut(1);
        while (FlagTimer2End == 0);
        FlagTimer2End = 0;
        DisableTimer2();
    }
}

void SimpleConfig(int time, int index) {
    initTimer2(time);
    EnableTimer2();
    QuickSet(index);
    while (FlagTimer2End == 0);
    FlagTimer2End = 0;
    DisableTimer2();
    
    int blank_state[11] = {0};

    initTimer2(1); //smallest period possible
    EnableTimer2();
    QuickSet(2);
    while (FlagTimer2End == 0);
    FlagTimer2End = 0;
    DisableTimer2();
}

void DirectConfig(int time, int sel, int val)
{
    initTimer2(10);
    EnableTimer2();
    DirectSet(sel, val);
    while (FlagTimer2End == 0);
    FlagTimer2End = 0;
    DisableTimer2();
    
    //int blank = 0b00000;

    initTimer2(1); //smallest period possible
    EnableTimer2();
    DirectSet(sel, 0); //blank, all switches off (shoot-through protection)
    while (FlagTimer2End == 0);
    FlagTimer2End = 0;
    DisableTimer2();
}

//Timer 1 interrupt 
void __attribute__((__interrupt__, no_auto_psv)) _T1Interrupt(void)
{
    //U1TXREG = 0x55; // Transmit one character
    FlagTimer1End = 1;
    IFS0bits.T1IF = 0; // clear TIMER5 interrupt flag
    //if (LATDbits.LATD4 == 0) LATDbits.LATD4 = 1;
    //else if (LATDbits.LATD4 == 1) LATDbits.LATD4 = 0;
}

//Timer 2 interrupt 
void __attribute__((__interrupt__, no_auto_psv)) _T2Interrupt(void)
{
    //U1TXREG = 0x55; // Transmit one character
    FlagTimer2End = 1;
    IFS0bits.T2IF = 0; // clear TIMER5 interrupt flag
    //if (LATDbits.LATD4 == 0) LATDbits.LATD4 = 1;
    //else if (LATDbits.LATD4 == 1) LATDbits.LATD4 = 0;
}
