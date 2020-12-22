#ifndef DATA_STRUCT_H_
#define DATA_STRUCT_H_


struct fluidCell {
   float ed, sd, temperature, pressure;
   float vx, vy, vz;
   float pi[4][4];
   float bulkPi;
};


struct fluidCell_2D {
    float temperature;
    float ux, uy, ueta;
    // the shear stress tensor are already divided by e+P
    float pi00, pi01, pi02;
    float pi11, pi12;
    float pi22;
    float pi33;
    float bulkPi;
};


struct fluidCell_3D_ideal {
    int itau, ix, iy, ieta;
    float energy_density;
    float temperature;
    float ux, uy, ueta;
};


struct fluidCell_3D {
    float temperature;
    float ux, uy, ueta;
    // the shear stress tensor are already divided by e+P
    float pi00, pi01, pi02, pi03;
    float pi11, pi12, pi13;
    float pi22, pi23;
    float pi33;
    float bulkPi;
};


struct fluidCell_3D_new {
    int itau, ix, iy, ieta;
    float temperature;
    float ux, uy, ueta;
    // the shear stress tensor are already divided by e+P
    float pi11, pi12, pi13;
    float pi22, pi23;
    float bulkPi;
};

#endif  // DATA_STRUCT_H_
