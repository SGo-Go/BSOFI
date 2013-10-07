#include <stdlib.h>
#include <string.h>

void system_(char * string, int len)
{ 
        char *str;

        str = malloc(len+1);
        strncpy(str,string,len);
        str[len]='\0';
        system(str);
}  
