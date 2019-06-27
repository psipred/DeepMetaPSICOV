/* Program to plot contact map using PSICOV output */

#include "gd.h"
#include "gdfontg.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include <string>
//  g++ plotmap.c -lgd -o plotmap
// gcc plotmap.c -lgd -o plotmap
using namespace std;

#define MAXSEQLEN 5000

void usage(){

	fprintf(stderr, "Usage: plotmap [options] <alnfile> <native contacts> <predicted contacts>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "-s	<float>	Score threshold. Default 0\n");
	fprintf(stderr, "-h	<int>	Show help\n\n");
}

int main(int argc, char **argv){
	//return(1);
	FILE *ifp;
	char seq[MAXSEQLEN],ss[MAXSEQLEN],line[MAXSEQLEN];
	int black,white,green,red,blue,r1,r2,d1,d2;
	float score,dist;
	float score_threshold = 0;
	int i = 1;
	int j = 0;
	int count = 0;

	// Buffer around grid
	int y_buffer = 160;
	int x_buffer = 160;
	// Size of each square
	int pixels_per_res = 30;
	// Text offsets
	int horizontal_text_x_offset = 10;
	int horizontal_text_y_offset = -25;
	int vertical_text_x_offset = -18;
	int vertical_text_y_offset = 7;
	int horizontal_text_x_scale_offset = 32;
	int horizontal_text_y_scale_offset = -25;
	int vertical_text_x_scale_offset = -30;
	int vertical_text_y_scale_offset = 32;

	int input_files[3];

	if (argc < 4) {
		usage();
		return 1;
	}

	while(i < argc){
                if(argv[i][0] == '-'){
                        i++;
                        switch(argv[i-1][1]){
				case 's' : {score_threshold = atof(argv[i]);
			printf("Score theshold is %f\n",score_threshold);break;}
                                case 'h' : {usage();return 1;}
                                default:   {usage();return 1;}
                        }
                }else{
			//printf("file %d: %s\n",i,argv[i]);
			input_files[j] = i;
			j++;
		}
		i++;
	}

	// Parse the alignment file to get sequence & sequence length
	ifp = fopen(argv[input_files[0]], "r");
	if (!ifp){
		fprintf(stderr, "Couldn't open alignment file\n");
		return 1;
	}

	fgets(seq, MAXSEQLEN, ifp);
	fclose(ifp);
	printf("Parsed alignment file.\n");

	int seq_len = strlen(seq)-1;
	//printf("Sequence is length %d\n",seq_len);

	// Some GD stuff
	// Declare the image
	gdImagePtr im;
	// Declare output files
	FILE *pngout, *jpegout;

	// Allocate the image
	im = gdImageCreate(x_buffer+seq_len*pixels_per_res,y_buffer+seq_len*pixels_per_res);

	// Allocate some colours
	white = gdImageColorAllocate(im, 255, 255, 255);
	black = gdImageColorAllocate(im, 0, 0, 0);
	green = gdImageColorAllocate(im, 5, 255, 20);
	red = gdImageColorAllocate(im, 220, 15, 15);
	blue = gdImageColorAllocate(im, 0, 0, 255);

	// Allocate some storage for contacts & scores
	int *n1_array = (int *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (int));
	if (n1_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}
	int *n2_array = (int *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (int));
	if (n2_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}


	// Parse native contacts
	ifp = fopen(argv[input_files[1]], "r");
	if (!ifp){
		fprintf(stderr, "Couldn't open native contacts\n");
		return 1;
	}

	// Parse observed contacts file
	while(fgets(line, 100, ifp)){
		if (sscanf(line, "%d%d%f", &r1, &r2, &dist)){
			if((r2-r1)>5){
				n1_array[count] = r1;
				n2_array[count] = r2;
				count++;
			}
		}
	}
	fclose(ifp);
	printf("Parsed native contacts.\n");

	// Draw contacts
//	for(i = 0; i < count; i++){
//		gdImageFilledRectangle(im,(x_buffer/2)+(n1_array[i]-1)*pixels_per_res,(y_buffer/2)+(n2_array[i]-1)*pixels_per_res,pixels_per_res+(x_buffer/2)+(n1_array[i]-1)*pixels_per_res,pixels_per_res+(y_buffer/2)+(n2_array[i]-1)*pixels_per_res,red);
//
//		gdImageFilledRectangle(im,(x_buffer/2)+(n2_array[i]-1)*pixels_per_res,(y_buffer/2)+(n1_array[i]-1)*pixels_per_res,pixels_per_res+(x_buffer/2)+(n2_array[i]-1)*pixels_per_res,pixels_per_res+(y_buffer/2)+(n1_array[i]-1)*pixels_per_res,red);
//	}
	printf("Drawn %d native contacts.\n",count);

	// Parse predicted contacts
	ifp = fopen(argv[input_files[2]], "r");
	if (!ifp){
		fprintf(stderr, "Couldn't open predicted contacts\n");
		return 1;
	}

	// Allocate some storage for contacts & scores
	int *r1_array = (int *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (int));
	if (r1_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}
	int *r2_array = (int *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (int));
	if (r2_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}

	// Allocate some storage for contacts & scores
	int *m1_array = (int *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (int));
	if (m1_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}
	int *m2_array = (int *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (int));
	if (m2_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}

	float *score_array = (float *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (float));
	if (score_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}

	int *ss1_array = (int *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (int));
	if (ss1_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}
	int *ss2_array = (int *) malloc((MAXSEQLEN*MAXSEQLEN) * sizeof (int));
	if (ss2_array == NULL) {
		printf("malloc failed\n");
		return 1;
	}

	char tmp[10];
	count = 0;
	// Parse predicted contacts file
	while(fgets(line, 100, ifp)){
		if (sscanf(line, "%d%d%d%d%f", &r1, &r2, &d1, &d2, &score)){
		//if (sscanf(line, "%d%d%d%f%[^,]", &r1, &r2, &d1, &score, tmp)){

			//printf("-%c--\n",tmp[1]);
			//if(tmp[1] == 'T'){
				//printf("----%c--\n",tmp[1]);
				if(score>=0.5){
					r1_array[count] = r1;
					r2_array[count] = r2;
					//score_array[count] = score;
					count++;
				}
			//}
		}
	}
	fclose(ifp);
	printf("Parsed predicted contacts.\n");

	// Draw contacts
	for(i = 0; i < count; i++){
		gdImageFilledRectangle(im,((x_buffer/2)+(r1_array[i]-1)*pixels_per_res)+2,((y_buffer/2)+(r2_array[i]-1)*pixels_per_res)+2,(pixels_per_res+(x_buffer/2)+(r1_array[i]-1)*pixels_per_res)-2,(pixels_per_res+(y_buffer/2)+(r2_array[i]-1)*pixels_per_res)-2,green);

		gdImageFilledRectangle(im,((x_buffer/2)+(r2_array[i]-1)*pixels_per_res)+2,((y_buffer/2)+(r1_array[i]-1)*pixels_per_res)+2,(pixels_per_res+(x_buffer/2)+(r2_array[i]-1)*pixels_per_res)-2,(pixels_per_res+(y_buffer/2)+(r1_array[i]-1)*pixels_per_res)-2,green);
	}
	printf("Drawn %d predicted contacts.\n",count);

	for(i = 0; i <= seq_len; i++){

		// Draw scale
		if(i % 25 == 0){

			// Draw hoizontal scale
			gdImageLine(im, (x_buffer/4), (y_buffer/2)+i*pixels_per_res, (x_buffer/4)+seq_len*pixels_per_res, (y_buffer/2)+i*pixels_per_res,red);

			// Draw vertical scale
			gdImageLine(im, (x_buffer/2)+i*pixels_per_res, (y_buffer/4), (x_buffer/2)+i*pixels_per_res,(y_buffer/4)+seq_len*pixels_per_res, red);

			char str[10];
			sprintf(str,"%d",i);

			gdImageString(im, gdFontGetGiant(), horizontal_text_x_scale_offset+(x_buffer/4)+i*pixels_per_res, horizontal_text_y_scale_offset+(y_buffer/4), (unsigned char*)str, black);

			gdImageString(im, gdFontGetGiant(), vertical_text_x_scale_offset+(x_buffer/4), vertical_text_y_scale_offset+(y_buffer/4)+i*pixels_per_res, (unsigned char*)str, black);

		}

		// Draw hoizontal lines
		gdImageLine(im, (x_buffer/2), (y_buffer/2)+i*pixels_per_res, (x_buffer/2)+seq_len*pixels_per_res, (y_buffer/2)+i*pixels_per_res, black);

		// Draw vertical lines
		gdImageLine(im, (x_buffer/2)+i*pixels_per_res, (y_buffer/2), (x_buffer/2)+i*pixels_per_res,(y_buffer/2)+seq_len*pixels_per_res, black);


		// Draw contacts for 'identical' residues
		if(i < seq_len){
			gdImageFilledRectangle(im,(x_buffer/2)+i*pixels_per_res,(y_buffer/2)+i*pixels_per_res,pixels_per_res+(x_buffer/2)+i*pixels_per_res,pixels_per_res+(y_buffer/2)+i*pixels_per_res,black);
		}

		// Draw text
		if(i < seq_len){

			gdImageChar(im, gdFontGetGiant(), horizontal_text_x_offset+(x_buffer/2)+i*pixels_per_res, horizontal_text_y_offset+(y_buffer/2), seq[i], black);

			gdImageChar(im, gdFontGetGiant(), vertical_text_x_offset+(x_buffer/2), vertical_text_y_offset+(y_buffer/2)+i*pixels_per_res, seq[i], black);

		}

	}

	// Open PNG file
	char* basename;
	char png[] = ".png";
	basename = strtok(argv[input_files[0]], ".");
	strcat(basename,png);
	printf("Writing output file %s...\n",basename);
	pngout = fopen(basename, "wb");

	// Write in PNG format
	gdImagePng(im, pngout);

	// Close PNG file
	fclose(pngout);

	// Free up memory
	gdImageDestroy(im);
	free(r1_array);
	r1_array = NULL;
	free(r2_array);
	r2_array = NULL;
	free(score_array);
	score_array = NULL;

	return 0;
}
