
typedef struct Spring {
    double k;
    float l0;
    int m1;
    int m2;
    float a;
    float b;
    float c;
} __attribute__ ((aligned (32))) Spring;

__kernel void breathing_cube(
    __global Spring* input,
    const float omega,
    const float t,
    const int data_size) {
    
    int i = get_global_id(0);
    if (i < data_size) {
        //printf("before in kernel: %f\n", input[i].l0);
        input[i].l0 = input[i].a + input[i].b * sin(omega*t + input[i].c);
        //printf("after in kernel: %f\n", input[i].l0);
    }
}

