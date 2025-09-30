/*
 * File:   main_PWMoutput.c
 * Author: Yiou He
 * Company: Massachusetts Institute of Technology, Power Electronics Research Group
 * Date: June 22, 2016
 * File Version:  0.1
 * Other Files Required: p33FJ64GS606.h
 * Tools Used: MPLAB ICD3
 * Devices Supported by this file: dsPIC 33FJ64GS606 
 */

//The HVPC can be started/stopped either by the RC controller or through the UART
//Through UART: send 'a' to precharge the input cap; send 'b' to start the HV; send 'c' to stop the HV
//All the parameters (desired output voltage, ADC calibration, etc) are set in main_PWMoutputDefinitions.h file
//This file records inflight data to FRAM through SPI;
//the "xxx_readclearSPI" program reads the data and save it in a txt file; 
//the "xxx_clearSPI" program clears the FRAM on board, prep for a new run.

/* File Description: */
// <editor-fold defaultstate="collapsed" desc="File Description">
/******************************************************************************* 
 * This program is for the high voltage power supply version 2 PCB board, it realizes the following functionalities:
 *  - Feedback control and PWM generation: 
 *    Measure the output voltage, input voltage and input current to adjust operating frequency, the controller goal is to fix the output voltage at 40 kV. 
 *  - Emergency shutdown when:
 *    Battery under voltage;
 *    Battery over temperature;
 *    RC control signal;
 *    Timer is up;
 *  - Storage data:
 *    Input voltage, current, output voltage data;
 *    Battery temperature data;
 *    Accelerometer data;
 * Requires to use the following function blocks:
 * - UART
 * - ADC
 * - PWM
 * - SPI                
********************************************************************************
 * Referenced files:
 * - Microchip example CE159-PWM_PushPull
 * - Microchip example CE160_StdPWM
 * - Microchip example CE166_Dual_Trig_ADC
 * - Dave Otten's assembly file ICNCmain_DaveOtten.s
 * - Microchip dsPIC33FJ64GS606 datasheet
 * - Microchip MPLABICD3 User's Guide for MPLAB X IDE          
********************************************************************************/
// </editor-fold>

#include "xc.h"
#include "p33FJ64GS606.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "dsPIC33FJ64GS606UART.h"
#include "dsPIC33FJ64GS606ADC.h"
#include "dsPIC33FJ64GS606SPI.h"
#include "dsPIC33FJ64GS606PWM.h"
#include "dsPIC33FJ64GS606CLK.h"
#include "dsPIC33FJ64GS606TIMER.h"
#include "dsPIC33FJ64GS606IO.h"
#include "main_PWMoutputDefinitions.h"
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
#pragma config HYST0 = HYST0
#pragma config HYST1 = HYST0
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

//Definitions
int CheckElectrical(double,double,double,double,double);  //Check to see if all values are OK
//int CheckMechamical(double,double,double,double);  //Check to see if all values are OK
void AdjustPWM(double,double,double,double);
//Operation functions
void Prepare(); //200V on, precharge capacitors
void Charge(); //charge
void Connect();  //200V on, connect converter with 200V
void HVOperation(); //200V on, HV on
void EndOperation(); //turn off everything
//define look up table
const double LookUpInputV[7] = {40.0,60.0,80.0,100.0,120.0,140.0,160.0};
const double LookUpOutputV[20] = {DesiredOutputV/20.0*10.0,DesiredOutputV/20.0*11.0,DesiredOutputV/20.0*12,
DesiredOutputV/20.0*13.0,DesiredOutputV/20.0*14,DesiredOutputV/20.0*14.5,DesiredOutputV/20.0*15.0,DesiredOutputV/20.0*15.4,DesiredOutputV/20.0*15.8,
DesiredOutputV/20.0*16.2,DesiredOutputV/20.0*16.6,DesiredOutputV/20.0*17.0,DesiredOutputV/20.0*17.4,DesiredOutputV/20.0*17.8,DesiredOutputV/20.0*18.2,
DesiredOutputV/20.0*18.6,DesiredOutputV/20.0*19.0,DesiredOutputV/20.0*19.4,DesiredOutputV/20.0*19.8,DesiredOutputV};

/*
const double LookUpOutputV[20] = {DesiredOutputV/20.0*17.5,DesiredOutputV/20.0*17.5,DesiredOutputV/20.0*17.5,
DesiredOutputV/20.0*18.0,DesiredOutputV/20.0*18.0,DesiredOutputV/20.0*18.0,DesiredOutputV/20.0*18.0,DesiredOutputV/20.0*18.5,DesiredOutputV/20.0*18.5,
DesiredOutputV/20.0*18.5,DesiredOutputV/20.0*18.5,DesiredOutputV/20.0*19.0,DesiredOutputV/20.0*19.0,DesiredOutputV/20.0*19.0,DesiredOutputV/20.0*19.0,
DesiredOutputV/20.0*19.5,DesiredOutputV/20.0*19.5,DesiredOutputV/20.0*19.5,DesiredOutputV/20.0*19.5,DesiredOutputV};
*/
const double LookUpFreq[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
int counter = 0;
const int maxcounter = 100;  //150 is 1.5 seconds
int IndexOutput = 0;
const int MaxIndexOutput = 19;
//Define PWM parameters
const int initPWMperiod = (int)((1000000/initPWMTimerFreq - 2000)/1.04) + 1880;
const int initPWMduty = (int)((((1000000/initPWMTimerFreq - 2000)/1.04) + 1880)/2.0);
int PWMperiod = (int)((1000000/initPWMTimerFreq - 2000)/1.04) + 1880; //normal operation
int MaxPWMperiod = (int)((1000000/MinPWMTimerFreq - 2000)/1.04) + 1880;
int MinPWMperiod = (int)((1000000/MaxPWMTimerFreq - 2000)/1.04) + 1880;
int PWMduty = (int)((((1000000/initPWMTimerFreq - 2000)/1.04) + 1880)/2.0); //normal operation
int MaxPWMduty = (int)((((1000000/MinPWMTimerFreq - 2000)/1.04) + 1880)/2.0);
int MinPWMduty = (int)((((1000000/MaxPWMTimerFreq - 2000)/1.04) + 1880)/2.0);
int PWMdeadtime = (int) initPWMdeadtime; //normal operation
//Define SPI and UART parameters
int long loc = 0;
unsigned int SPIdata = 0;
int SPIdataCounter = 0;
int SPIdataRecorder = 1;
char buf[30];
//Define ADC parameters
int ADCperiod = ADCtimerPeriod*Timer12Freq;  //ADC sample and conversion once in ADCperiod/460.5 ms
int ADCcounter = 0;
int ADCsampletime = 0;
double ADCFloat[5] = {0.0,0.0,0.0,0.0,0.0};
double InputVoltage = 0;
double InputCurrent = 0;
double Temperature = 0;
double OutputVoltage = 0;
double LogicBV = 0;
double AccX = 0;
double AccY = 0;
double AccZ = 0;
//Define flag for error code
int FlagE = 0; //Electrical value check error code
//int FlagM = 0; //Mechanical value check error code
int IC1cnt = 0;
unsigned int IC1timeAve = 0;
unsigned int IC1timeTotal = 0;
int FlagINCON = 0; //input capture interrupt switch ON
//int FlagINCOFF = 0; //input capture interrupt switch OFF
int FlagTimerEnd = 0; //timer end, then ignore RC control
//Define timer end interrupt
int EndTimeMSB = 0xD2;  //End time in sec, 0xD2 for 29.49 MHz (115.1953*1000*120), 0x1A5 is for 29.49*2 MHz (230.3906*1000*120 = 27646872)
int EndTimeLSB = 0xCA80;  //End time in sec, 0xCA80 for 29.49 MHz, 0xDB98 for 29.49*2 MHz

int main()
{
    //Initialization
    initIO();
    DischargeCaps();
    initCLK();
    initPWM(initPWMperiod, 0, initPWMdeadtime); //initiate with 0% duty cycle for PWM1, 100% for PWM2
    initUART1();
    //initUART2();
    initSPI2();
    //initial ADC
    initTimer1(ADCperiod);
    initADC(initPWMperiod);
    //initial timer counter for 2 mins
    initTimer45(EndTimeMSB,EndTimeLSB);
    //initial input interrupt for RC signal
    initTimer2();
    initIC1();
    EnableTimer1();
    EnableTimer2();
    EnableTimer45();
    //EnableADC();
    EnablePWM();
    EnableSPI2();
    EnableUART1();
    WriteScreen("HelloWorld");
    sprintf(buf,"MSB is %d , LSB is %u ",EndTimeMSB,EndTimeLSB);
    WriteScreen(buf);
    WriteScreen("\n\r");
    while(1)
    {
        if (ADCsampletime == 10)
        {
            //LATFbits.LATF6 = 1;
            InputVoltage = ADCFloat[0];
            InputCurrent = ADCFloat[3];///ADCsampletime;
            OutputVoltage = ADCFloat[2];///ADCsampletime;
            Temperature = ADCFloat[1];///ADCsampletime;
            LogicBV = ADCFloat[4];///ADCsampletime;
            ADCsampletime = 0;
            
            FlagE = CheckElectrical(InputVoltage,InputCurrent,OutputVoltage,LogicBV,Temperature);  //Check to see if all values are OK
            //Noticed that when start up, the input current is always big, so can we bypass it?
            if (FlagE == 0)
            {
                if (counter == maxcounter)
                {
                    counter = 0;
                    IndexOutput = IndexOutput + 1;
                }
                else counter = counter + 1;
                if (IndexOutput >= MaxIndexOutput) IndexOutput = MaxIndexOutput;
                AdjustPWM(InputVoltage,InputCurrent,OutputVoltage,LookUpOutputV[IndexOutput]);
            }
            else 
            {
                EndOperation();
                DisableTimer45();
            }
            if (SPIdataCounter == 50)  //every 10*50*ADCtiming, SPI
            {
                SPIdataCounter = 0;
                SPIdataRecorder ++;
            }
            loc = SendThroughSPI(SPIdataRecorder,InputVoltage,InputCurrent,OutputVoltage,LogicBV,loc);
            //WriteScreenADCdata(FlagE,InputVoltage,InputCurrent,OutputVoltage,Temperature,LogicBV,AccX,AccY,AccZ);
            SPIdataCounter ++;
            
            //Reset ADC
            ADCFloat[0]=0.0;
            ADCFloat[1]=0.0;
            ADCFloat[2]=0.0;
            ADCFloat[3]=0.0;
            ADCFloat[4]=0.0;
            //LATFbits.LATF6 = 0;
        }
    }
}

// Capture Interrupt Service Routine
void __attribute__((__interrupt__, no_auto_psv)) _IC1Interrupt(void)
{
    unsigned int t1,t2;
    unsigned int IC1timePeriod;
    //LATDbits.LATD4 = 1;
    t1 = IC1BUF;
    t2 = IC1BUF;
    IFS0bits.IC1IF=0;
    if(t2 > t1)
        IC1timePeriod = t2 - t1;
    else
        IC1timePeriod = (PR2 - t1) + t2;
    if (IC1timePeriod <= HVOFFminperiod) return;
    if (IC1timePeriod >= HVONmaxperiod) return;
    //Take Average
    if (IC1cnt == 5)
    {
        IC1timeAve = floor(((float)(IC1timeTotal))/5.0);
        IC1timeTotal = 0;
        IC1cnt = 0;
    }
    else
    {
        IC1cnt = IC1cnt + 1;
        IC1timeTotal = IC1timeTotal + IC1timePeriod;
    }
    //calculate see whether is turn-on signal or turn-off signal
    if (IC1timeAve >= HVONminperiod)
    {
        if (IC1timeAve <= HVONmaxperiod)
        {
            if (FlagTimerEnd == 0)
            {
            switch(FlagINCON)
            {
                case 0:
                    Prepare();
                    FlagINCON = 1;
                    //WriteScreen("Pre-charge Caps\n\r");
                    //PreparePWM();
                    //Connect();
                    //HVOperation();
                    break;
                case 1:
                    FlagINCON = FlagINCON + 1;
                    //PreparePWM();
                    break;
                case 10:  //50
                    FlagINCON = FlagINCON + 1;
                    //Connect();
                    break;
                case 60:
                    FlagINCON = FlagINCON + 1;
                    PreparePWM();
                    Connect();
                    HVOperation();
                    break;
                default:
                    FlagINCON = FlagINCON + 1;
            }
            }
        }
    }
    //Following section is correct, does not cause current overshoot
    if (IC1timeAve >= HVOFFminperiod)
    {
        if (IC1timeAve <= HVOFFmaxperiod)
        {
            FlagINCON = 0;
            FlagTimerEnd = 0;
            EndOperation();
            DisableTimer45();
        }
    }
    if (FlagTimerEnd == 1)
    {
        FlagINCON = 0;
        EndOperation();
        DisableTimer45();
    }
}

//Timer 5 interrupt 
void __attribute__((__interrupt__, no_auto_psv)) _T5Interrupt(void)
{
    //U1TXREG = 0x55; // Transmit one character
    FlagTimerEnd = 1;
    EndOperation(); //timer is up
    //DisableTimer45();
    IFS1bits.T5IF = 0; // clear TIMER5 interrupt flag
    //if (LATDbits.LATD4 == 0) LATDbits.LATD4 = 1;
    //else if (LATDbits.LATD4 == 1) LATDbits.LATD4 = 0;
}

void __attribute__((__interrupt__, no_auto_psv)) _U1RXInterrupt(void) {

    char receive = 0;
    receive = U1RXREG; // Receive one character
    //U1TXREG = receive;
    
    switch (receive)
    {
        case 'a':
            U1TXREG = 'a'; //Letter 'a'
            WriteScreen("Pre-charge Caps\n\r");
            Prepare();
            break;
        case 'b':
            U1TXREG = 'b';
            WriteScreen("Normal Operation\n\r");
            PreparePWM();
            Connect();
            HVOperation();
            break;
        case 'c':
            EndOperation();
            U1TXREG = 'c';
            WriteScreen("Converter stopped\n\r");
            break;
        case 'd':
            U1TXREG = 'd';
            WriteScreen("Converter running\n\r");
            break;
        case '=':
            PWMperiod = PWMperiod + 2;
            PWMduty = PWMduty + 1;
            ChangePD(PWMperiod, PWMduty, initPWMdeadtime);
            WriteScreen("Increase period\n\r");
            break;
        case '-':
            PWMperiod = PWMperiod - 2;
            PWMduty = PWMduty - 1;
            ChangePD(PWMperiod, PWMduty, initPWMdeadtime);
            WriteScreen("Decrease period\n\r");
            break;
        case 'w':
            PWMdeadtime = PWMdeadtime + 2;
            ChangePD(PWMperiod, PWMduty, PWMdeadtime);
            WriteScreen("Increase deadtime\n\r");
            break;
        case 'e':
            PWMdeadtime = PWMdeadtime - 2;
            ChangePD(PWMperiod, PWMduty, PWMdeadtime);
            WriteScreen("Decrease deadtime\n\r");
            break;
        default:
            break;
    }
    char buf[50];
    sprintf(buf,"Period is %d, duty is %d, dtime is %d    ",PWMperiod, PWMduty, PWMdeadtime);
    WriteScreen(buf);
    WriteScreen("\n\r");
    IFS0bits.U1RXIF = 0; // clear RX interrupt flag
}

//UART transmitter interrupt, worked. 
void __attribute__((__interrupt__, no_auto_psv)) _U1TXInterrupt(void)
{
    IFS0bits.U1TXIF = 0; // clear TX interrupt flag
    //U1TXREG = 0x55; // Transmit one character
}

//UART transmitter interrupt, worked.
void __attribute__((__interrupt__, no_auto_psv)) _U2TXInterrupt(void)
{
    IFS1bits.U2TXIF = 0; // clear TX interrupt flag
    //U2TXREG = 0x55; // Transmit one character
}

//ADC interrupt

void __attribute__((__interrupt__, no_auto_psv)) _ADCP0Interrupt(void) {
    IFS6bits.ADCP0IF = 0;
    //LATFbits.LATF6 = 1;
    ADCFloat[0] = (double)((unsigned int)(7.0*ADCFloat[0] + (double)ADCBUF0)>>3);//((double) ADCBUF0) * (double) (3.3/1023);
    ADCFloat[1] = (double)((unsigned int)(7.0*ADCFloat[1] + (double)ADCBUF1)>>3);//((double) ADCBUF1) * (double) (3.3/1023);
    //LATFbits.LATF6 = 0;
}

void __attribute__((__interrupt__, no_auto_psv)) _ADCP1Interrupt(void) {
    IFS6bits.ADCP1IF = 0;
    ADCFloat[2] = (double)((unsigned int)(7.0*ADCFloat[2] + (double)ADCBUF2)>>3);
    ADCFloat[3] = (double)((unsigned int)(7.0*ADCFloat[3] + (double)ADCBUF3)>>3);
}

void __attribute__((__interrupt__, no_auto_psv)) _ADCP2Interrupt(void) {
    IFS7bits.ADCP2IF = 0;
    ADCFloat[4] = (double)((unsigned int)(7.0*ADCFloat[4] + (double)ADCBUF4)>>3);
    //ADCFloat[5] = (double) ADCBUF5 + ADCFloat[5];
    ADCsampletime = ADCsampletime + 1;
}
/*void __attribute__((__interrupt__, no_auto_psv)) _ADCP3Interrupt(void) {
    IFS7bits.ADCP3IF = 0;
    ADCFloat[6] = (double) ADCBUF6 + ADCFloat[6];
    ADCFloat[7] = (double) ADCBUF7 + ADCFloat[7];
    //LATFbits.LATF6 = 0;
    ADCsampletime = ADCsampletime + 1;
}
*/
//SPI transmit interrupt
void __attribute__((__interrupt__, no_auto_psv)) _SPI2Interrupt(void)
{
    //U2TXREG = SPIdata; // Transmit the received SPIdata
    //if (LATFbits.LATF6 == 1)
    //    LATFbits.LATF6 = 0;
    //else LATFbits.LATF6 = 1;
    IFS2bits.SPI2IF = 0; // clear SPI2 interrupt flag
}

int CheckElectrical(double InputVoltage,double InputCurrent,double OutputVoltage,double LogicBV, double Temperature)  //Check to see if all values are OK
{
    if (IndexOutput == 0)
        return 0;
    else
    {
        if (((InputVoltage*InputVScale + InputVOffset) < MinInputV) || ((InputVoltage*InputVScale + InputVOffset) > MaxInputV)) return 1;
        else if ((InputCurrent*InputIScale + InputIOffset) > MaxInputI) return 2;
        //else if ((OutputVoltage*OutputVScale + OutputVOffset) > MaxOutputV) return 3;
        //else if (((LogicBV*LogicBVScale + LogicBVOffset) < MinLogicBV) || ((LogicBV*LogicBVScale + LogicBVOffset) > MaxLogicBV)) return 4;
        //else if (Temperature*TempScale + TempOffset > MaxTemp) return 5;
        return 0;
    }
}

void AdjustPWM(double InputVoltage,double InputCurrent,double OutputVoltage,double SetOutputV)
{
    if ((OutputVoltage*OutputVScale + OutputVOffset) < SetOutputV)
    {
        PWMperiod = PWMperiod + 2;
        PWMduty = PWMduty + 1;
    }
    else if ((OutputVoltage*OutputVScale + OutputVOffset) > SetOutputV)
    {
        PWMperiod = PWMperiod - 2;
        PWMduty = PWMduty - 1;
    }
    if (PWMperiod > MaxPWMperiod) PWMperiod = MaxPWMperiod;
    if (PWMperiod < MinPWMperiod) PWMperiod = MinPWMperiod;
    if (PWMduty > MaxPWMduty) PWMduty = MaxPWMduty;
    if (PWMduty < MinPWMduty) PWMduty = MinPWMduty;
    ChangePD(PWMperiod, PWMduty, initPWMdeadtime);
    //Check with look up table
}

void Prepare() //200V on, precharge capacitors
{
    //clear timer45 flag
    FlagTimerEnd = 0;
    DisablePWM();
    DisableTimer45();
    //Change PWM to initial states
    initChangePD(initPWMperiod, 0, initPWMdeadtime);
    PWMperiod = initPWMperiod;
    PWMduty = initPWMduty;
    PWMdeadtime = initPWMdeadtime;
    //restart output voltage profile and create delays
    counter = 0;
    IndexOutput = 0;
    //PrechargePWM();
    //PreparePWM();
    Precharge();
}

void Connect()  //200V on, connect converter with 200V
{
    EnablePWM();
    //Precharge();
}
void HVOperation() //200V on, HV on
{
    ChangePD(initPWMperiod, initPWMduty, initPWMdeadtime);
    NormalOperation();
    HVon();
    //EnablePWM();
    EnableTimer45(); //start timer 45 to detect 2 mins.
    EnableADC();
    EnableSPI2();
}
void EndOperation()
{
    //clear timer flag
    //FlagTimerEnd = 0;
    DischargeCaps();
    //Reinitiate PWM
    DisablePWM();
    //PrechargePWM();
    DisableADC();
    DisableSPI2();
    counter = 0;
    IndexOutput = 0;
} 