/* 
 * File:   proj2.c
 * Author: Luká¹ Vokráèko
 *
 * Created on October 16, 2012, 10:15 PM
 * 
 * Description:
 *	Výpoèet a^x, tan^-1(x), argsinh(x) se zadanou pøesnotí sigdig
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

//prototypy fcí
void printMsg(const char** msgArr, int msgType);
void printResult(double num);
void readInput(const int function, const double a, const double epsilon);

bool getNumParam(char* param, double* ptrnum, bool readint);
bool isPrecisionValid(double sum, double num, const double epsilon);
bool isAbout(double num, int val, const double epsilon);

int getParam(char* param);

double customAtan(double x, const double epsilon);
double _customAtan(double x, const double epsilon);
double customArgsinh(double x, const double epsilon);
double customLn(double x, const double epsilon);
double _customLnE(double x, const double epssylon);
double customPow(double base, double exp, const double epsilon);
double customExp(double base, const double epsilon);
double getEpsilonFromSigdig(const int sigdig);

//definice èísel s nekoneèným rozvojem
const double IZPLN2 = 0.6931471805599453094172321214581765680755001343602552;
const double IZPPIPUL = 1.5707963267948966192313216916397514420985846996875529;

//definice indexù parametrù
enum PARAM_TYPES
{
    EXP,
    ARCTG,
    ARGSINH,
    HELP,
    NOPARAM,
};

//definice indexù chyb
enum ERROR_TYPES
{
    EPARAM,
    ESIGDIGHIGH,
};

//definice parametrù
const char* PARAMS[] ={
    [EXP] = "--powxa",
    [ARCTG] = "--arctg",
    [ARGSINH] = "--argsinh",
    [HELP] = "-h",
};

//definice chyb
const char* ERRORS[] = {
    [EPARAM] = "Chybnì zadané parametry\nPro nápovìdu spus»te program s parametrem -h",
    [ESIGDIGHIGH] = "Se zadanou pøesností není mo¾né pracovat, nastavuji na 15",
};

//definice stringu obsahujícího nápovìdu
const char* HELPSTR =
	"NÁZEV\n"
	"	proj2\n"
	"AUTOR\n"
	"	Luká¹ Vokráèko\n"
	"POU®ITÍ\n"
	"	proj2 [--powxa|--argsinh|--arctg|-h] sigidg a\n"
	"PARAMETRY\n"
	"	-h - tiskne nápovìdu\n"
	"	--powxa sigdig a - vypoèet x^a\n"
	"	--arctg sigidg - výpoèet tan^-1(x)\n"
	"	--argsinh sigdig - výpoèet argsinh(x)\n"
	"	sigdig [1-15] poèet platných èíslic\n"
	"	a - exponent mocninné funkce\n"
	;

//zji¹tìní zadaného parametru
int getParam(char* param)
{
    if(strcmp(param, PARAMS[EXP]) == 0) return EXP;
    else if(strcmp(param, PARAMS[ARCTG]) == 0) return ARCTG;
    else if(strcmp(param, PARAMS[ARGSINH]) == 0) return ARGSINH;
    else if(strcmp(param, PARAMS[HELP]) == 0) return HELP;
    else return NOPARAM;
}

//pøeètení èísla ze stringu - buï double nebo int podle parametru readint
bool getNumParam(char* pargvpos, double* ptrnum, bool readint)
{
    char* endptr;
    //ètu double
    if(readint == false) 
    {
	*ptrnum = strtod(pargvpos, &endptr);
    }
    //ètu int a pøevádím do double
    else
    {
	*ptrnum = (double) strtol(pargvpos, &endptr, 10);
    }
    if(*endptr == '\0')
    {
	return true;
    }
    return false;
}

//je hodnota témìø pøesná?
bool isAbout(double num, int val, const double epsilon)
{
    return (num >= val-epsilon && num <= val+epsilon); 
}

//fce pro poèítání exponencionální fce
double customPow(double base, double exp, const double epsilon)
{
    double integral_part;
    double integral_res = 1;
    
    double fractional_part;
    double fractional_res = 1;
    double res;
    
    //o¹etøení definièního oboru
    if(exp == 1) return base;
    if(exp == -1) return 1/base;
    if(base == 1 || exp == 0) return 1;
    if(isinf(base) && isinf(exp)) return INFINITY;
    if(isnan(base) || isnan(exp)) return NAN;
    
    if(base == 0)
    {
	if(exp == 0) return NAN;
	else if(exp > 0 || isinf(exp) == 1) return 0;
	else return INFINITY;
    }
    
    if(exp == 0)
    {
	if(isinf(base) != 0) NAN;
	else return 1;
    }
    
    if(isinf(exp) == -1) 
    {
	if(base == -1) return 1;
	else return 0;
    }
    if(isinf(exp) == 1)
    {
	if(base == -1) return 1;
	else return INFINITY;
    }
    
    if(isinf(base) != 0)
    {
	if(exp < 0) return 0;
	else return INFINITY;
    }
    
    //rozdìlení exponentu na celou a desetinnou èást
    fractional_part = modf(fabs(exp), &integral_part);
    if(base < 0 && fractional_part != 0) return NAN;
    
    //výpoèet celé èásti exponentu
    for(int i = 0; i < integral_part; i++)
    {
	integral_res *= base;
    }
    
    //výpoèet desetinné èásti exponentu
    if(fractional_part != 0) 
    {
	fractional_part *= customLn(base, epsilon);
	fractional_res = customExp(fractional_part, epsilon);
    }
    
    res = integral_res * fractional_res;
    
    //záporný exponent
    if(exp < 0) return 1/res;
    return res;
}


//fce pro výpoèet e^x
double customExp(double x, const double epsilon)
{
    bool neg  = (x < 0) ? true : false;
    x = fabs(x);
    double next = x;
    int i = 1;
    double sum = 1 + next;
    
    //výpoèet øady
    while(!isPrecisionValid(sum, next, epsilon))
    {
	++i;
	next = next * (x/i);
	sum += next;
    }
    
    if(neg == true) return 1/sum;
    else return sum;
}

//je pøesnost platná?
bool isPrecisionValid(double sum, double num, const double epsilon)
{
     return (fabs(num) <= fabs(sum*epsilon));
}

//vlastní fce pro výpoèet pøirozeného logaritmu
double customLn(double x, const double epsilon)
{
    double res;
    int twoTimes;
    
    //pro extrémní hodnoty není tøeba pokou¹et se zmìn¹ovat
    if(isinf(x)) return INFINITY;
    else if(x == 0) return -INFINITY;
    else if(x < 0 || isnan(x)) return NAN;

    //dìlení dvìmi dokud se nedostaneme do intervalu <0.5;1)
    x = frexp(x, &twoTimes);
    //o¹etøení definièního obru
    
    if(isAbout(x, 1, epsilon)) res = 0;
    //výpoèet
    else
    {
	res = _customLnE(x, epsilon);
    }
    //matematické úpravy pro získání správného výsledku
    res += twoTimes*IZPLN2;
    
    return res;
}

//výpoèet pøirozeného logartimu
double _customLnE(double x, const double epsilon)
{
    int k = 2;
    double s = (-1)*(-1+x);
    double sumup = s*s;
    double next = sumup/k;
    double sum = s + next;

    //výpoèet øady
    while(!isPrecisionValid(sum, next, epsilon))
    {
	++k;
	sumup *= s;
	next = sumup/k;
	sum += next; 
    }
    
    return -sum;
}

//zji¹tìní kolikrát lze èíslo vydìlit dvìmi a modifikace èísla
int getDivisionsByTwo(double* x)
{
    int counter = 0;
    while(*x >= 2)
    {
	*x /= 2;
	++counter;
    }
    while(*x < 1 )
    {
	*x *= 2;
	--counter;
    }
    return counter;
}

//vlastní fce pro výpoèet arctgx
double customAtan(double x, const double epsilon)
{
    if(x == 0) return 0;
    if(isinf(x) != 0) 
    {
	if(isinf(x) == 1) return IZPPIPUL;
	else return -IZPPIPUL;
    }
    if(isnan(x)) return NAN;
    
    double res = 0;
    double minusone = 1;
    
    if(x < 0) 
    {
	minusone = -1;
	x = fabs(x);
    }
    //matematická úprava
    if(x > IZPPIPUL)
    {
	res = IZPPIPUL - _customAtan(1/x, epsilon);
    }
    else
    {
	res = _customAtan(x, epsilon);
    }
    return res*minusone;
    
}

//fce pro výpoèet arctgx v intervalu (0;pi/2)
double _customAtan(double x, const double epsilon)
{
    double xx =x*x;
    double y = xx/(1+xx);
    double k = 2;
    double prev = 1;
    double next = (2*y)/3;
    double sum = prev + next;
    
    //výpoèet øady
    while(!isPrecisionValid(sum, next, epsilon))
    {
	k += 2;
	next *= ((k*y)/(k+1));
	sum += next;
    }
    
    return (y/x)*sum;    
}

//fce pro výpoèet argsinh
double customArgsinh(double x, const double epsilon)
{
    if(isinf(x) == -1) return -INFINITY;
    double res = x+customPow(x*x+1, 0.5, epsilon);
    return customLn(res, epsilon);
}

//formát pro výpis výsledku
void printResult(double num)
{
    printf("%.10e\n", num);
}

//výpis chybových hlá¹ení
void printMsg(const char** msgArr, int msgType)
{
    fprintf(stderr, "%s\n", msgArr[msgType]);
}

//pøevedení sigdig na epsilon
double getEpsilonFromSigdig(const int sigdig)
{
    return customPow(10, -sigdig, 1);
}

//ètení posloupnosti èísel na stdin a provádìní zadaných fcí
void readInput(const int function, const double a, const double epsilon)
{
    double num;

    int code;
    while((code = scanf("%lf", &num)) != EOF || code == 0)
    {
	//¹patný vstup, bereme jako nan
	if(code == 0) {
	    scanf("%*s");
	    num = NAN;
	}
	
	//jakou fci pou¾ít
	switch(function)
	{
	    case EXP:
		num = customPow(num, a, epsilon);
		break;
	    case ARGSINH:
		num = customArgsinh(num, epsilon);
		break;
	    case ARCTG:
		num = customAtan(num, epsilon);
		break;
	}
	
	//tisk výsledku
	printResult(num);
    }
}


int main(int argc, char** argv)
{   
    int param;
    double sigdig;
    double a;
    
    //naètení parametru
    if(argc >= 2)
    {
	 param = getParam(argv[1]);
    }
    
    //param -h?
    if((argc == 2 && param == HELP))
    {
	printf("%s\n", HELPSTR);
    }
    //jiný parametr, pøeètení pøesnosti
    else if((argc == 3 || argc == 4) && getNumParam(argv[2], &sigdig, true))
    {
	if(sigdig <= 0)
	{
	    printMsg(ERRORS, EPARAM);
	    return (EXIT_FAILURE);
	}
	else if(sigdig > DBL_DIG)
	{
	    printMsg(ERRORS, ESIGDIGHIGH);
	    sigdig = DBL_DIG;
	}
	const double epsilon = getEpsilonFromSigdig(sigdig+1);
	
	if(param == ARCTG && argc == 3)
	{
	    readInput(ARCTG, 0, epsilon);
	}
	else if(param == ARGSINH && argc == 3)
	{
	    readInput(ARGSINH, 0, epsilon);
	}
	//--powxa má 1 parametr navíc
	else if(param == EXP && argc == 4 && getNumParam(argv[3], &a, false))
	{
	    readInput(EXP, a, epsilon);
	}
	else 
	{
	    //neznámý parametr / ¹patný parametr pro a
	    printMsg(ERRORS, EPARAM);
	    return (EXIT_FAILURE);	    
	}
    }
    //¹patný poèet parametrù / ¹patný parametr pro sigdig
    else
    {
	printMsg(ERRORS, EPARAM);
	return (EXIT_FAILURE);
    }
    
    return (EXIT_SUCCESS);
}

