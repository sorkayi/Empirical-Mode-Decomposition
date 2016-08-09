/*
 * The MIT License
 *
 * Copyright (c) 2010 Hian-Kun Tenn
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * -----------------------------------------------------------------------------
 * This code is used to test the EMD functions
 * written by kylecimcdonald.
 *
 * The original project is held at:
 * http://code.google.com/p/realtime-emd/
 *
 * Author:      Hian-Kun Tenn
 * Date:        2010/12/29
 */
#include <stdio.h>
#include "EmpiricalModeDecomposition.h"

void imfs_output(emdData* emd){
        int i, j;

        for(j = 0; j < emd->size; j++){
                for(i = 0; i < emd->order; i++) {
                        printf("%f ", emd->imfs[i][j]);
                }
                printf("\n");
        }
}

int main(void){
        FILE *inputfile = fopen("sinusoid.txt","r");
        int i=0;
        float signal[1024]={0.0};
        emdData data;

        while (fscanf(inputfile, "%f", &signal[i]) != EOF){
                //printf("%d %f\n",i,signal[i]);
                i++;
        }

        emdCreate(&data,i,5,10,0);
        emdDecompose(&data,signal);
        imfs_output(&data);
        emdClear(&data);

        fclose(inputfile);
        return 0;
}
