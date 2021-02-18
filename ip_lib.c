/*ID Gruppo 45 E banalmente vero Danci Romina Alina 881597,Basile Giorgio 882779,Saran Mattia 881900,Perencin Francesco 880106*/ 
/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

#define E 2.7182818284590

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i < a->h && j < a->w &&k < a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}

/**** CREATE ****/

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
    ip_mat *a;
    int i,j,z;
    
    
    a = (ip_mat *) malloc (sizeof(ip_mat));
    
    if(a == NULL){
        printf("Errore ip_mat_create");
        exit(1);
    }
    
    a->h = h;
    a->w = w;
    a->k = k;
    
    a->data = (float ***) malloc (sizeof(float**) * h);
    if(a->data == NULL){
        printf("Errore ip_mat_Create");
        exit(1);
    }
    
    for(i=0;i<h;i++){
        a->data[i] = (float **) malloc (sizeof(float *) * w);
        if(a->data[i] == NULL){
            printf("Errore ip_mat_Create");
            exit(1);
        }
        for(j=0;j<w;j++){
            a->data[i][j] = (float *) malloc (sizeof(float) * k);
                if(a->data[i][j] == NULL){
                    printf("Errore ip_mat_Create");
                    exit(1);
                }
        }
   }

   for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            for(z=0;z<k;z++)
                set_val(a,i,j,z,v);
            
   a->stat = (stats *) malloc (sizeof(stats) * k);
    if(a->stat == NULL){
        printf("Errore ip_mat_Create");
        exit(1);
    }
   
   for(i=0;i<k;i++){
       a->stat[i].min = 0;
       a->stat[i].max = 0;
       a->stat[i].mean = 0;
   }
   
return a;
}

/**** FREE ****/

void ip_mat_free(ip_mat *a){
    int i,j;
    
    if(a != NULL){
        for(i=0;i < a->h;i++){
            for(j=0;j < a->w;j++){
                
                free(a -> data[i][j]);
                
            }
        }
        
        
        for(i=0;i < a->h;i++)
            free(a->data[i]);
        
        if(a->data != NULL)
            free(a->data);
        
        if(a->stat != NULL)
            free(a->stat);
        
        free(a);
    }
}

/**** COMPUTE  ****/

void compute_stats(ip_mat * t){
    int i,j;
    float somma;
    float celle;
    float min;
    float max;
    int can = 0;
    
    if(t == NULL){
        printf("Errore ip_mat_compute");
        exit(1);
    }
    
    while(can < 3){
        min = (get_val(t,0,0,can));
        max = (get_val(t,0,0,can));
        
        celle = 0.0;
        somma = 0.0;
        
        for(i=1;i<t->h;i++){
            for(j=1;j<t->w;j++){
                somma = somma + (get_val(t,i,j,can));
                if((get_val(t,i,j,can)) > max)
                    max = (get_val(t,i,j,can));
                if((get_val(t,i,j,can)) < min)
                    min = (get_val(t,i,j,can));
                celle = celle + 1;
            }
        }
        
        t->stat[can].min = min;
        t->stat[can].max = max;
        t->stat[can].mean = somma / celle;
        
        can = can + 1;
    }
                
}

/**** INNIT RANDOM ****/

void ip_mat_init_random(ip_mat * t, float mean, float var){
    
    int i,j,z;

    if(t == NULL){
        printf("Errore ip_mat_innit_random");
        exit(1);
    }
    
    for(i=0; i < t->h; i++){
         for(j=0; j < t->w; j++){
             for(z=0; z < t->k; z++){
                 
                 float m;
                 
                 m = get_normal_random(mean, var);
                 set_val(t,i,j,z,m);
             }
         }
    }
}

/**** COPY ****/

ip_mat * ip_mat_copy(ip_mat * in){
    ip_mat * out;
    float v;
    int i,j,z;

    if(in == NULL){
        printf("Errore ip_mat_copy");
        exit(1);
    }
    
    out = ip_mat_create(in->h, in->w, in->k, 0.0);
    
   for(i=0;i < in->h;i++){
        for(j=0;j < in->w;j++){
            for(z=0;z < in->k;z++){
                
                v = get_val(in,i,j,z);
                set_val(out,i,j,z,v);

            }
        }
   }
            
   return out;
}

/**** SUBSET ****/

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){

    ip_mat* a;
    int i,j,z;
    float m;
    
    if(t == NULL){
        printf("Errore ip_mat_subset");
        exit(1);
    }
    
    a = ip_mat_create(row_end - row_start, col_end - col_start, t->k, 0.0);

    for(i = row_start; i < row_end; i++){
         for(j = col_start; j < col_end; j++){
             for(z = 0; z < t -> k; z++){
                 
                 m = get_val(t,i,j,z);
                 set_val(a,i-row_start,j-col_start,z,m);
                 
             }
         }
    }

    return a;
}

/**** CONCAT ****/

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
    
    ip_mat*out;
    int i,j,z;
    float m;
    
    if(a == NULL || b == NULL || dimensione > 2){
        printf("Errore ip_mat_concat");
        exit(1);
    }
    
    if(dimensione == 0){
        
        out = ip_mat_create(a->h+b->h,a->w,a->k,0.0);
        
        for(i=0 ; i<out->h ; i++){
            for(j=0 ; j<out->w ; j++){
                for(z=0 ; z<out->k ; z++){
                    if(i<a->h){
                        m=get_val(a,i,j,z);
                        set_val(out,i,j,z,m);
                    } else {
                         m=get_val(b,i-(a->h),j,z);
                         set_val(out,i,j,z,m);
                    }
                }
            }
        }
    } else if(dimensione == 1){
        
        out=ip_mat_create(a->h,a->w+b->w,a->k,0.0);
        
        for(i=0 ; i<out->h ; i++){
            for(j=0 ; j<out->w ; j++){
                for(z=0 ; z<out->k ; z++){
                    if(j<a->w){
                        m=get_val(a,i,j,z);
                        set_val(out,i,j,z,m);
                    } else {
                         m=get_val(b,i,j-(a->w),z);
                         set_val(out,i,j,z,m);
                    }
                }
            }
        }
    } else if(dimensione == 2){
        
        out=ip_mat_create(a->h,a->w,a->k+b->k,0.0);
        
        for(i=0 ; i<out->h ; i++){
            for(j=0 ; j<out->w ; j++){
                for(z=0 ; z<out->k ; z++){
                    if(z<a->k){
                        m=get_val(a,i,j,z);
                        set_val(out,i,j,z,m);
                    } else {
                         m=get_val(b,i,j,z-(a->k));
                         set_val(out,i,j,z,m);
                    }
                }
            }
        }
    } else {
        printf("Dati inseriti non correttamente");
        exit(1);
    }
    return out;
}

/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/

/*** SUM ***/

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    
    ip_mat* c;
    int i,j,z;
    float sum;
    
    if(a->w == b->w && a->h == b->h && a->k == b->k && a != NULL && b!= NULL){
        c = ip_mat_create(a->h,a->w,a->k,0.0);
    } 
    else {
        printf("Errore ip_mat_sum");
        exit(1);
    }
    
    for(i=0 ; i < a->h ; i++){
        for(j=0 ; j < a->w ; j++){
            for(z=0 ; z < a->k ; z++){
                
                sum = get_val(a, i, j, z) + get_val(b, i, j, z);
                set_val(c, i, j, z, sum);
            }
        }
    }
    
    return c;
}

/*** SUB ***/

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    
    ip_mat* c;
    int i,j,z;
    float sub;
    
    if(a->w == b->w && a->h == b->h && a->k == b->k && a != NULL && b != NULL){
        c = ip_mat_create(a->h,a->w,a->k,0.0);
    } 
    else{
        printf("Errore ip_mat_sub");
        exit(1);
    }
    for(i=0 ; i < a->h ; i++){
        for(j=0 ; j < a->w ; j++){
            for(z=0 ; z < a->k ; z++){
                
                sub = get_val(a, i, j, z) - get_val(b, i, j, z);
                set_val(c, i, j, z, sub);

            }
        }
    }
    
    return c;
}


/*** MUL SCALAR ***/

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    
    int i,j,z;
    float mul;
    ip_mat* d;
    
    if(a == NULL){
        printf("Errore ip_mat_mul_scalar");
        exit(1);
    }
    
    d = ip_mat_create(a->h, a->w, a->k, 0.0);
    
    for(i=0 ; i < a->h ; i++){
        for(j=0 ; j < a->w ; j++){
            for(z=0 ; z < a->k ; z++){
                
                mul = get_val(a, i, j, z) * c;
                set_val(d, i, j, z, mul);
                
            }
        }
    }
    
    return d;
}

/*** ADD SCALAR ***/

ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
    
    int i,j,z;
    float sum;
    ip_mat* d;
    
    if(a == NULL){
        printf("Errore ip_mat_add_scalar");
        exit(1);
    }
    
    d = ip_mat_create(a->h, a->w, a->k, 0.0);
    
    for(i=0 ; i < a->h ; i++){
        for(j=0 ; j < a->w ; j++){
            for(z=0 ; z < a->k ; z++){
                
                sum = get_val(a, i, j, z) + c;
                set_val(d, i, j, z, sum);
                
            }
        }
    }
    
    return d;
}

/*** MEAN ***/

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){

    ip_mat* c;
    
    if(a->w == b->w && a->h == b->h && a->k == b->k && a != NULL && b != NULL){
        c = ip_mat_create(a->h, a->w, a->k, 0.0);
    } 
    else {
        printf("Errore ip_mat_mean");
        exit(1);
    }
    
    c = ip_mat_sum(a,b);
    c = ip_mat_mul_scalar(c,0.5);
    
    return c;
}


/************ 2 PARTE ************/

/***TO GRAY SCALE***/

ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    
    ip_mat* out = ip_mat_create(in->h,in->w,in->k,0.0);
    int i,j,z;
    float mean;
    
    if(in == NULL){
        printf("Errore ip_mat_gray_scale");
        exit(1);
    }
    
    for(i=0 ; i<out->h ; i++){
        for(j=0 ; j<out->w ; j++){
            mean = (get_val(in, i, j, 0) + get_val(in, i, j, 1) + get_val(in, i, j, 2)) / 3;
            for(z=0; z<out->k; z++){
                
                set_val(out, i, j, z, mean);
                
            }
        }
    }
    
    return out;
}

/***BRIGHTEN***/

ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    
    ip_mat* result;
    
    if(a == NULL){
        printf("Errore ip_mat_brighten");
        exit(1);
    }
    
    result = ip_mat_add_scalar(a, bright);
    
    return result;
}

/***BLEND***/

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    ip_mat* primo;
    ip_mat* secondo;
    ip_mat* blend;
    
    
    if(a == NULL && b == NULL){
        printf("Errore ip_mat_blend");
        exit(1);
    }
    
    primo = ip_mat_mul_scalar(a,alpha);
    secondo = ip_mat_mul_scalar(b,(1.0-alpha));
    blend = ip_mat_sum(primo,secondo);
    ip_mat_free(primo);
    ip_mat_free(secondo);
    
    return blend;
    
}

/***CORRUPT***/

ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    
    ip_mat* t;
    
    if(a == NULL){
        printf("Errore ip_mat_corrupt");
        exit(1);
    }

    t = ip_mat_create(a->h,a->w,a->k,0.0);
    
    ip_mat_init_random(t, 0, 1.);
    t = ip_mat_mul_scalar(t, amount);

    t = ip_mat_sum(a,t);
    
    return t;
}

/************ 3 PARTE ************/

/**FILTRO SHARPEN**/

ip_mat * create_sharpen_filter(){
    ip_mat *s;
    int i,j;
    int r = 0;
    float val[9] = {0.0,-1.0,0.0,-1.0,5.0,-1.0,0.0,-1.0,0.0};
    
    s = ip_mat_create(3, 3, 1, 0.0);
    
    if(s == NULL){
        printf("Errore ip_mat_sharpen_filter");
        exit(1);
    }
    
    for(i = 0;i < s->h; i++){
        for(j = 0;j < s->w; j++){
            set_val(s,i,j,0,val[r]);
            r++;
        }
    }
    
    return s;
}

/**FILTRO EDGE**/

ip_mat * create_edge_filter(){
    ip_mat *ed;
    int i,j;
    int r = 0;
    float val[9] = {-1.0,-1.0,-1.0,-1.0,8.0,-1.0,-1.0,-1.0,-1.0};
    
    ed = ip_mat_create(3, 3, 1, 0.0);
    
    if(ed == NULL){
        printf("Errore ip_mat_edge_filter");
        exit(1);
    }
    
    for(i = 0;i < ed->h; i++){
        for(j = 0;j < ed->w; j++){
            set_val(ed,i,j,0,val[r]);
            r++;
        }
    }
    
    return ed;
}

/**FILTRO EMBOSS**/

ip_mat * create_emboss_filter(){
    ip_mat *em;
    int i,j;
    int r = 0;
    float val[9] = {-2.0,-1.0,0.0,-1.0,1.0,1.0,0.0,1.0,2.0};
    
    em = ip_mat_create(3, 3, 1, 0.0);
    
    if(em == NULL){
        printf("Errore ip_mat_emboss_filter");
        exit(1);
    }
    
    for(i = 0;i < em->h; i++){
        for(j = 0;j < em->w; j++){
            set_val(em,i,j,0,val[r]);
            r++;
        }
    }
    
    return em;
}

/**FILTRO GAUSSIAN**/

ip_mat * create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma){
    ip_mat * out;
    int x, y;
    int i,j,z;
    float somma = 0.0;
    int cx, cy;
    float esp, val,base;
    
    out = ip_mat_create(h, w, k, 0.);
    
    if(out == NULL){
        printf("Errore ip_mat_gaussian_filter");
        exit(1);
    }
    
    cx = h / 2;
    cy = w / 2;
    
    for(i = 0;i < h; i++){
        for(j = 0;j < w; j++){
            for(z = 0;z < k;z++){
                
                x = i - cx;
                y = j - cy;
                
                esp = -(((pow(x,2) + (pow(y,2))) / (2*(pow(sigma,2)))));
                base = (1 / (2 * PI * (pow(sigma,2))));
                val = pow(base,esp);
                
                set_val(out,i,j,z,val);
                
                somma = somma + val;
            }
        }
    }
    
    out = ip_mat_mul_scalar(out, (1 / somma));
    
    return out;
}

/**FILTRO AVERAGE**/

ip_mat * create_average_filter(unsigned w, unsigned h, unsigned k){ 
    ip_mat *av;
    
    float c = 1.0 /((float) (w*h));
    
    av = ip_mat_create(w, h, k, c);
    
    if(av == NULL){
        printf("Errore ip_mat_averege_filter");
        exit(1);
    }
    
    return av;
}

/**PADDING**/

ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w){
    
    int i,j,z;
    float m;
    ip_mat * out;
    
    out = ip_mat_create(a->h+2*pad_h,a->w+2*pad_w,a->k,0.0);
    
    if(a == NULL || out == NULL){
        printf("Errore ip_mat_padding");
        exit(1);
    }
    
    for(i = pad_h;i < ((out->h) - pad_h); i++){
        for(j = pad_w; j < ((out->w)-pad_w); j++){
            for(z = 0; z < out->k;z++){
                
                m = get_val(a, i-pad_h, j-pad_w, z);
                set_val(out, i, j, z, m);
                
            }
        }
    }
    
    return out;
}

/**RESCALE**/

void rescale(ip_mat * t, float new_max){
    int i,j,can;
    float r;
    
    if(t == NULL){
        printf("Errore ip_mat_rescale");
        exit(1);
    }
    
    compute_stats(t);
    
    for(can=0; can < 3; can++){
        for(i=0;i < t->h;i++){
            for(j=0;j < t->w;j++){
                r = get_val(t, i, j, can) - (t->stat[can].min) / ((t->stat[can].max) - (t->stat[can].min));
                set_val(t, i, j, can, r);
            }
        }
    }
    
    t = ip_mat_mul_scalar(t, new_max);
}

/**CLAMP**/

void clamp(ip_mat * t, float low, float high){
    
    int i,j,z;
    
    if(t == NULL){
        printf("Errore ip_mat_clamp");
        exit(1);
    }
    
    rescale(t,255.0);
    
    for(i=0;i<t->h;i++){
        for(j=0;j<t->w;j++){
            for(z=0;z<t->k;z++){

                if(get_val(t, i, j, z) < low)
                    set_val(t, i, j, z, low);
                else{
                    if(get_val(t, i, j, z) > high)
                        set_val(t, i, j, z, high);
                }
            }
        }
    }
}

/**CONVOLVE**/

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    int i, j,z;
    int ph, pw;
    ip_mat * out;
    ip_mat * b;
    int gauss = 0;
        
    if(a == NULL || f == NULL){
        printf("Errore ip_mat_convolve");
        exit(1);
    }

    ph = ((f->h)-1)/2;
    pw = ((f->w)-1)/2;

    out = ip_mat_padding(a, ph, pw);
    b = ip_mat_copy(a);
    
    if(out == NULL || b == NULL){
        printf("Errore ip_mat_mul_convolve");
        exit(1);
    }
    
    if((get_val(f,0,0,0)) != (get_val(f,1,1,0))){
        gauss = 1;
    }

    for(i = 0;i < b->h; i++){
        
        for(j = 0;j < b->w; j++){
            
            int i1,j1,z1;
            ip_mat * t;
            t =  ip_mat_subset(out, i, i + (f->h), j, j + (f->w));
            
            for(z=0;z<3;z++){
                
                float somma=0.0;
                
                if(f->k == 3 && gauss){
                        for(i1 = 0;i1 < f->h; i1++){
                            for(j1 = 0;j1 < f->w; j1++){
                                for(z1 = 0;z1 < f->k;z1++){
                                    somma = somma + ((t->data[i1][j1][z])*(f->data[i1][j1][z1]));
                                }
                            }
                        }
                        set_val(b, i, j, z, somma);
                }
                else{
                    for(i1 = 0;i1 < f->h; i1++){
                        for(j1 = 0;j1 < f->w; j1++){
                            somma = somma + ((t->data[i1][j1][z])*(f->data[i1][j1][0]));
                        }
                    }
                    set_val(b, i, j, z, somma);
                }
            }
            ip_mat_free(t);
        }
    }
    
    return b;
}


ip_mat * bitmap_to_ip_mat(Bitmap * img){
    
    unsigned int i=0,j=0;
    
    unsigned char R;
    unsigned char G;
    unsigned char B;
    
     
    int h = img->h;
    int w = img->w;
    
    
    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    compute_stats(out);

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}
