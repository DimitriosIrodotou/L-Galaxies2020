# 1 "/Users/Bam/ClionProjects/Rings/h_metals.h"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 361 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "/Users/Bam/ClionProjects/Rings/h_metals.h" 2

/*struct elements_str
{
    float H;
    float He;
    float O;
    float Mg;
    float Fe;
};
#endif  //MAINELEMENTS*/
/* Define two views into the same element structure/array.
   For example define:
      union elements DiskMassElements;
   Then access as either
      DiskMassElements.str.H
   or
      DiskMassElements.arr[0]
*/
/*
union elements
{
    struct elements_str str;
    float arr[NUM_ELEMENTS];
};*/



//Number of chemical elements tracked:
