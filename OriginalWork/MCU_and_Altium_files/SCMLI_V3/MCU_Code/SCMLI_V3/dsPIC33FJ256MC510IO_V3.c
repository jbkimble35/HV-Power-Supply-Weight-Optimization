/*
 * File:   dsPIC33FJ64GS606IO.c
 * Author: Suzanne
 *
 * Created on July 7, 2016, 9:53 AM
 */

#include "dsPIC33FJ256MC510IO_V3.h"

//bits from R to L: SWE' to SWA
//states: 0 = charge, 1 = +0, 2 = +1, 3 = -1
int Gen1 = 0b11001;
int Gen2 = 0b11000;
int Gen3 = 0b01010;
int Gen4 = 0b10100;
int Top1 = 0b00000;
int Top2 = 0b10000;
int Top3 = 0b00010;
int Top4 = 0b10000;
int Bot1 = 0b01001;
int Bot2 = 0b01000;
int Bot3 = 0b01000;
int Bot4 = 0b00100;
int All0 = 0b00000;
int Boot = 0b01010; 

//which letter down, which module across
//mod1 is A2,C1,D1
//mod2 is A3,B2,C2,D2,E1
int port_indices[5][11] = {
    {3,5,0,6,4,2,0,1,0,5,7}, //A2 to A12
    {7,3,6,4,4,4,6,4,1,5,5}, //B1 to B11
    {0,3,5,4,4,4,1,2,4,1,7}, //C1 to C11
    {0,3,0,6,6,2,6,1,1,1,7}, //D1 to D11
    {7,3,3,6,6,4,2,6,1,1,1}, //E0 to E10
};

//which letter down, which module across
int bit_indices[5][11] = {
    {1, 0, 7,13,5, 3,0,3,1,3, 17}, //A2 to A12
    {16,3, 1,1, 3, 7,6,9,1,12,2},   //B1 to B11
    {14,2, 1,0, 2, 6,6,4,8,2, 17},   //C1 to C11
    {15,13,6,12,15,2,8,4,6,13,17},   //D1 to D11
    {17,8,12,0, 14,4,7,5,1,0, 12}    //E1 to E11
};

//config down, module across
int *sw_configs[6][11] = {
    {&All0,&All0,&All0,&All0,&All0,&All0,&All0,&All0,&All0,&All0,&All0},
    {&Bot1,&Gen1,&Gen1,&Gen1,&Gen1,&Gen1,&Gen1,&Gen1,&Gen1,&Gen1,&Top1},
    {&Bot2,&Gen2,&Gen2,&Gen2,&Gen2,&Gen2,&Gen2,&Gen2,&Gen2,&Gen2,&Top2},
    {&Bot3,&Gen3,&Gen3,&Gen3,&Gen3,&Gen3,&Gen3,&Gen3,&Gen3,&Gen3,&Top3},
    {&Bot4,&Gen4,&Gen4,&Gen4,&Gen4,&Gen4,&Gen4,&Gen4,&Gen4,&Gen4,&Top4},
    {&Boot,&Boot,&Boot,&Boot,&Boot,&Boot,&Boot,&Boot,&Boot,&Boot,&Boot}
};

int bit_configs[20][7] = {{0}};

//placeholder variables used to assign values to the ports
int port_configs[7] = {
    0x0000,
    0x0000,
    0x0000,
    0x0000,
    0x0000,
    0x0000,
    0x0000
};

void initIO()  //working, tested on 7/4/2016
{
    //Initialize the IO port
    //Use PORT to read, use LAT to write
    //PORTA
    LATA = 0x0000;  //RD: many SW outputs, 1 is INT_3
    TRISA = 0x0000; //0b0000000000000000;
    //PortB 
    LATB = 0x0000;  //Set the latch to zeros for all
    TRISB = 0x0000; //0b0000000000000000;
    //PORTC 
    LATC = 0x0000;
    TRISC = 0x0000; //0b0000000000000000;
    //PORTD 
    LATD = 0x0000;  //RD: many SW outputs, 1 is INT_3
    TRISD = 0x0000; //0b0000000000000000;
    //PORTE 
    LATE = 0x0000;  //RE: many SW outputs, 5 is INT_5
    TRISE = 0x0000; //0b0000000000000000;
    //PORTF 
    LATF = 0x0000;
    TRISF = 0x0000; //0b0000000000000000
    //PORTG 
    LATG = 0x0000;  //RG: 9 is INT, others are SW or nothing
    TRISG = 0x0000; //0b0000000000000000
    //Analog and comparator input setting
    //All analog and comparator pins are default as analog input, need to be set to 1 when it's used as digital pins
    //Not using any analog inputs as of now
    ADPCFG = 0xFFFF;
    //Disable comparator pairs
}

void SetModules(int *state_indices, int index) 
{
    int config;
    int port;
    int bito;
    int curr_por;
    int x;
    int y;
    for(x=0; x<11; x++)
    {
        config = *(sw_configs[*state_indices][x]);
        
        for(y=0; y<5; y++)
        {
            port = port_indices[y][x];
            bito = bit_indices[y][x];
            curr_por = port_configs[port];
            //figure out if you need to write a pin high for some config
            if (((config >> y) & 1U) && (port < 7) && (bito < 17))
            {
                curr_por |= 1UL << bito;
            }
            else
            {
                curr_por &= ~(1UL << bito);
            }
            port_configs[port] = curr_por;
        }
        /*increment pointer for next element fetch*/
        *state_indices++;
    }
    
    for (x=0; x<7; x++)
    {
        //port_configs[x] = ~port_configs[x];
        bit_configs[index][x] = port_configs[x];
    }
    WriteScreen(bit_configs);
}

void TestOut(int sel)
{
    if (sel == 0) {
        LATA = 0x0000;
        LATB = 0x0000;
        LATC = 0x0000;
        LATD = 0x0000;
        LATE = 0x0000;
        LATF = 0x0000;
        LATG = 0x0000;
    }
    else if (sel == 1) {
        LATA = 0xFFFF;
        LATB = 0xFFFF;
        LATC = 0xFFFF;
        LATD = 0xFFFF;
        LATE = 0xFFFF;
        LATF = 0xFFFF;
        LATG = 0xFFFF;
    }
    else if (sel == 2) {
        LATA = 0x5555;
        LATB = 0x5555;
        LATC = 0x5555;
        LATD = 0x5555;
        LATE = 0x5555;
        LATF = 0x5555;
        LATG = 0x5555;
    }
    else {
        LATA = 0xAAAA;
        LATB = 0xAAAA;
        LATC = 0xAAAA;
        LATD = 0xAAAA;
        LATE = 0xAAAA;
        LATF = 0xAAAA;
        LATG = 0xAAAA;
    }
}

void DirectSet(int sel, int val)
{
    //flip to save ur dumb butt
    //val = ~val;
        
    switch(sel) {
        case 0  :
           LATA = val;
           break; /* optional */
        case 1  :
           LATB = val;
           break; /* optional */
        case 2  :
           LATC = val;
           break; /* optional */
        case 3  :
           LATD = val;
           break; /* optional */
        case 4  :
           LATE = val;
           break; /* optional */
        case 5  :
           LATF = val;
           break; /* optional */
        case 6  :
           LATG = val;
           break; /* optional */
        /* you can have any number of case statements */
        default : /* Optional */
           LATA = val;
    }
}

void QuickSet(int index) 
{
    LATA = bit_configs[index][0];
    LATB = bit_configs[index][1];
    LATC = bit_configs[index][2];
    LATD = bit_configs[index][3];
    LATE = bit_configs[index][4];
    LATF = bit_configs[index][5];
    LATG = bit_configs[index][6];
}

int TestRead(int sel_int) {
    return PORTB;
}
