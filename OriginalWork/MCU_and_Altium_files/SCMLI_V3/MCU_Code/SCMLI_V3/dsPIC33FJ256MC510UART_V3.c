/*
 * File:   dsPIC33FJ64GS606UART.c
 * Author: Yiou
 *
 * Created on July 7, 2016, 9:34 AM
 */

#include "dsPIC33FJ256MC510UART_V3.h"

void initUART1()
{
    U1MODEbits.STSEL = 0; // 1 Stop bit
    U1MODEbits.PDSEL = 0; // No Parity, 8 data bits
    U1MODEbits.ABAUD = 0; // Auto-Baud Disabled
    U1MODEbits.BRGH = 0; // Low Speed mode
    
    U1BRG = 1;
    /*
     *  ;29.49 MHz Fcy (64/16 x PLL with 7.37MHz FRC)
		;191 = 9600 (+0.00%)
		;95 = 19200 (+0.00%)
		;47 = 38400 (+0.00%)
		;31 = 56000 (+0.00%)
		;15 = 115200 (+0.00%)
		;7 = 230400 (+0.00)
		;3 = 460800 (+0.00%)
		;1 = 921600 (+0.00%)
    */
    /*
     *  ;58.98 MHz Fcy (64/8 x PLL with 7.37MHz FRC
		;191 = 19200 (+0.00%)
		;95 = 38400 (+0.00%)
		;47 = 56000 (+0.00%)
		;31 = 115200 (+0.00%)
		;15 = 230400 (+0.00)
		;7 = 460800 (+0.00%)
		;3 = 921600 (+0.00%)
    */    
    U1STA = 0; //init status and control register, interrupt after one TX Character is transmitted
    U1TXREG = 0x0000;
    //Bit 15-9: Unimplemented (0000000)
    //Bit 8: Data bit 8 in 9-bit mode (0)
    //Bit 7-0: Data bits 7-0 (00000000)
       
    //Set interrupt
    IFS0bits.U1RXIF = 0;
    IFS0bits.U1TXIF = 0;
    IEC0bits.U1RXIE = 1;
    IEC0bits.U1TXIE = 1;
    //IPC2bits.U1RXIP = 0b111;
}
void EnableUART1()
{
    U1MODEbits.UARTEN = 1; //Enable UART
    Nop();
    U1STAbits.UTXEN = 1;   //Enable UART transmit
    Nop();
}

void DisableUART1()
{
    U1MODEbits.UARTEN = 0; //Enable UART
    Nop();
    U1STAbits.UTXEN = 0;   //Enable UART transmit
    Nop();
}

//UART2 for the bluetooth
void initUART2()
{
    U2MODEbits.STSEL = 0; // 1 Stop bit
    U2MODEbits.PDSEL = 0; // No Parity, 8 data bits
    U2MODEbits.ABAUD = 0; // Auto-Baud Disabled
    U2MODEbits.BRGH = 0; // Low Speed mode
    
    U2BRG = 15;
    /*
     *  ;29.49 MHz Fcy (64/16 x PLL with 7.37MHz FRC)
		;191 = 9600 (+0.00%)
		;95 = 19200 (+0.00%)
		;47 = 38400 (+0.00%)
		;31 = 56000 (+0.00%)
		;15 = 115200 (+0.00%)
		;7 = 230400 (+0.00)
		;3 = 460800 (+0.00%)
		;1 = 921600 (+0.00%)
    */
    
    U2STA = 0; //init status and control register, interrupt after one TX Character is transmitted
    U2TXREG = 0x0000;
    //Bit 15-9: Unimplemented (0000000)
    //Bit 8: Data bit 8 in 9-bit mode (0)
    //Bit 7-0: Data bits 7-0 (00000000)
       
    //Set interrupt
    IFS1bits.U2RXIF = 0;
    IFS1bits.U2TXIF = 0;
    IEC1bits.U2RXIE = 0;
    IEC1bits.U2TXIE = 1;
    //IPC2bits.U1RXIP = 0b111;
}

void EnableUART2()
{
    U2MODEbits.UARTEN = 1; //Enable UART
    Nop();
    U2STAbits.UTXEN = 1;   //Enable UART transmit
    Nop();
    LATBbits.LATB8 = 1;    //Set DSR pin high
}
void WriteScreen(char s[50]) 
{
    char *p;
    p = s;
    while (*p)
    {
        while (!(U1STAbits.TRMT)) {}
        U1TXREG = *(p++);
    }
}
void WriteBluetooth(char s[50]) 
{
    char *p;
    p = s;
    while (*p)
    {
        while (!(U1STAbits.TRMT)) {}
        U2TXREG = *(p++);
    }
}

void WriteScreenADCdata(int FlagE, double InputVoltage,double InputCurrent,double OutputVoltage,double Temperature,double LogicBV,double AccX,double AccY,double AccZ)
{
    char buf[3];
    sprintf(buf,"%d",FlagE);
    WriteScreen(buf);
    WriteScreen("\r\n");
    sprintf(buf,"%.2f",InputVoltage); //*InputVScale + InputVOffset
    WriteScreen(buf);
    WriteScreen("\r\n");
    sprintf(buf,"%.2f",InputCurrent); //*InputIScale + InputIOffset
    WriteScreen(buf);
    WriteScreen("\r\n");        
    sprintf(buf,"%.2f",OutputVoltage); //*OutputVScale + OutputVOffset
    WriteScreen(buf);
    WriteScreen("\r\n");
    sprintf(buf,"%.2f",LogicBV);//LogicBV*LogicBVScale + LogicBVOffset);
    WriteScreen(buf);
    WriteScreen("\r\n");
    sprintf(buf,"%.2f",Temperature);//LogicBV*LogicBVScale + LogicBVOffset);
    WriteScreen(buf);
    WriteScreen("\r\n");
}