//**************************************************************************************
//							THRAIN I/O
//  ThrainIO.cpp                                                         S Reece Boston
//	This is a first step towads an I/O interface.
//	Will take a filename form the commandline, uses as input
//		The input files need to be formatted in the expected way to function
//      See documentation, or sample input, to see how to use
//  Rece Boston, Mar 24, 2022
//**************************************************************************************

#include "ThrainMain.h"
#include "ThrainIO.h"

//this function reads in user input from the specified file to create calculation data
int read_input(char input_file_name[128], CalculationInputData &calcdata){
	//open file
	FILE* input_file;
	try {
			input_file = fopen(input_file_name, "r");
	} catch (std::exception& e) {
			printf("Input file not found.\n");
			return 1;
	};
	
	//The following properties of the calculation are defined in the input file
	char input_buffer[128];		//for read-in from file and text processing
	std::string instring;		//for more read-in from file
	char calculation_name[128];	//a label for this calculation
	regime::Regime regime;		//the physics to be used (Newtonian, 1PN, GR)
	model::StellarModel model;	//the model of star to be used
	units::Units units;			//the units of the calculation
	modetype::ModeType modetype;//the type of modes (radial, nonradial, cowling, etc.)
	calcdata.mass = calcdata.radius = calcdata.zsurf = calcdata.logg = calcdata.teff = 0.0;
	calcdata.params = 0;

	//We want to allow for comments to be added to the top of an input file
	//These will be preceded by a "#" symbol
	// NOTE: comments can ONLY go at the top of an input file
	int startofline=ftell(input_file);
	int lines=0;
	char chr = getc(input_file);
	printf("%c\n", chr);
	while(chr=='#' | chr=='\n'){
		if(chr=='#') fgets(input_buffer, 128, input_file);
		chr = getc(input_file);
		lines++;
	}
	fseek(input_file, startofline, SEEK_SET);
	for(int line=0; line<lines; line++) fgets(input_buffer, 128, input_file);
	
	
	//Now we shall extract info from the input file
	fscanf(input_file, "Name: %s\n", calculation_name);
	calcdata.calcname = std::string(calculation_name);
	printf("calculation name: %s\n", calcdata.calcname.c_str());
	
	//read the type of physics to be used in calculation
	fscanf(input_file, "Model: %s\t", input_buffer);
	instring = std::string(input_buffer);
	if(!instring.compare("newtonian")) calcdata.regime = regime::PN0;
	else{
		printf("ERROR IN INPUT: THRAIN only knows Newtonian physics!\n");
		printf("This error is recoverable!  Setting to proper regime.\n");
		calcdata.regime = regime::PN0;
	}
	//read the background stellar model to be used in calculation
	fscanf(input_file, "%s\t", input_buffer);
	printf(" MAKING A %s\n", input_buffer);
	instring = std::string(input_buffer);
	if(!instring.compare("polytrope")) calcdata.model = model::polytrope;
	else if(!instring.compare("CHWD")) calcdata.model = model::CHWD;
	else if(!instring.compare("MESA")) calcdata.model = model::MESA;
	else if(!instring.compare("SWD"))  calcdata.model = model::SWD;
	else {
		printf("ERROR IN INPUT: THRAIN can only work with polytrope, CHWD, MESA, or SWD\n");
		printf("This error is fatal.  Quitting.");
		return 1;
	}
	
	//handle the use of different parameters that help construct the stellar model
	//POLYTROPE INPUT
	if(calcdata.model==model::polytrope){
		calcdata.num_input_params = 2;
		calcdata.input_params = new double[calcdata.num_input_params];
		fscanf(input_file, " %lf", &calcdata.input_params[0]);	//read in the index
		printf("n=%lf\n", calcdata.input_params[0]);
		fscanf(input_file, "%lf\n", &calcdata.input_params[1]);	//read in the grid size
		calcdata.Ngrid = int(calcdata.input_params[1]);
		
		//now read in desired physical properties of star -- specify two at a time
		//for a polytrope, can specify two of mass, radius, logg, or surface z
		double temp;
		//read in first physical parameter -- name as string, value as double
		fscanf(input_file, "Params: %s %lf\t", input_buffer, &temp);
		instring=std::string(input_buffer);
		//save value in appropriate slot use bit-masking to keep track of variables
		if(     !instring.compare("mass"))   {calcdata.mass = temp; calcdata.params|=units::pmass;}
		else if(!instring.compare("radius")) {calcdata.radius=temp; calcdata.params|=units::pradius;}
		else if(!instring.compare("zsurf"))  {calcdata.zsurf =temp; calcdata.params|=units::pzsurf;}
		else if(!instring.compare("logg"))   {calcdata.logg  =temp; calcdata.params|=units::plogg;}
		char tempparam = calcdata.params;	//value after first read-in, for comparison later
		//read in second physical parameters -- name as string, value as double
		fscanf(input_file, "%s %lf\n", input_buffer, &temp);
		instring=std::string(input_buffer);
		//save valeu in appropriate slot, again use bit-masking so that binary value of params indicates which were used
		if(     !instring.compare("mass"))   {calcdata.mass = temp; calcdata.params|=units::pmass;}
		else if(!instring.compare("radius")) {calcdata.radius=temp; calcdata.params|=units::pradius;}
		else if(!instring.compare("zsurf"))  {calcdata.zsurf =temp; calcdata.params|=units::pzsurf;}
		else if(!instring.compare("logg"))   {calcdata.logg  =temp; calcdata.params|=units::plogg;}
		else {
			printf("ERROR IN INPUT: invalid parameter to polytope:\n allowed mass, radius, zsurf, or logg\n");
			printf("This error is fatal.  Quitting.\n");
			return 1;
		}
		//make sure we did not specify the same parameter twice
		if(calcdata.params == tempparam){
			printf("ERROR IN INPUT: double specification of stellar parameters, star underdefined\n");
			printf("This error is fatal.  Quitting.\n");
			switch(calcdata.params){
				case units::pmass:   printf("\tmass given twice\n"); break;
				case units::pradius: printf("\tradius given twice\n"); break;
				case units::pzsurf:  printf("\tzsurf given twice\n"); break;
				case units::plogg:   printf("\tlog g given twice\n"); break;
			}
			return 1;	
		}	
		if((calcdata.params&units::pteff)){
			printf("ERROR IN INPUT: polytropes cannot be assigned temperature\n");
			printf("This error is fatal.  Quitting.\n");
			return 1;
		}
	}
	
	//CHANDRASEKHAR WD INPUT
	else if(calcdata.model==model::CHWD){
		calcdata.num_input_params = 3;
		calcdata.input_params = new double[calcdata.num_input_params];
		fscanf(input_file, " %lf", &calcdata.input_params[0]);	//read in the central y value
		printf("y0=%lf\n", calcdata.input_params[0]);
		//read in an integer designating the chemical composition to use for mu
		fscanf(input_file, " %lf", &calcdata.input_params[1]);
		switch((int) calcdata.input_params[1]){
			case 0:
				printf("standard Chandrasekhar WD\n");
				break;
			case 1:
				printf("using logistic composition\n");
				break;
		}
		fscanf(input_file, "%lf\n", &calcdata.input_params[2]);	//read in the grid size
		calcdata.Ngrid = int(calcdata.input_params[2]);
		//fscanf(input_file, "%s\n", input_buffer);   //read in chemical gradient
		//calcdata.str_input_param = std::string(input_buffer);
	}
	
	//MESA INPUT
	else if(calcdata.model==model::MESA){
		calcdata.num_input_params = 1;
		calcdata.input_params = new double[calcdata.num_input_params];
		fscanf(input_file, "%s ", input_buffer);    //read in the filename
		calcdata.str_input_param=std::string(input_buffer);
		printf("source=%s.dat\n", calcdata.str_input_param.c_str());
		fscanf(input_file, "%lf\n", &calcdata.input_params[0]);	//read in the grid size
		calcdata.Ngrid = int(calcdata.input_params[0]);
	}
	
	//SIMPLE WD INPUT
	else if(calcdata.model==model::SWD){
		calcdata.num_input_params = 10;
		calcdata.input_params = new double[calcdata.num_input_params];
		
		fscanf(input_file, "%lf\n", &calcdata.input_params[0]);	//read in the grid size
		calcdata.Ngrid = int(calcdata.input_params[0]);
		printf("Ngrid=%d\n", calcdata.Ngrid);
		
		//on a new line, specify the EOS
		FILE *swd = fopen("swd.txt", "w");
		fscanf(input_file, "EOS:\n");
		fprintf(swd, "# equation of state\n");
		//read the EOS specification
		fscanf(input_file, "\tcore\t%[^\n]", input_buffer);
		fprintf(swd, "core:\n\t%s\n", input_buffer);
		calcdata.str_input_param = std::string(input_buffer);
		calcdata.str_input_param += std::string("\n");
		fscanf(input_file, "\tatm\t%[^\n]", input_buffer);
		fprintf(swd, "atm:\n\t%s\n\n", input_buffer);
		calcdata.str_input_param += std::string(input_buffer);
		
		//on a new line, specify chemical paramters
		fscanf(input_file, "chemical parameters:\n");  //skip a line of formatting
		fprintf(swd, "# chemical parameters\n");
		fscanf(input_file, " %*[^0123456789] %lf %lf %lf\n", &calcdata.input_params[1], &calcdata.input_params[2], &calcdata.input_params[3]);	//read in the helium
		fscanf(input_file, " %*[^0123456789] %lf %lf %lf\n", &calcdata.input_params[4], &calcdata.input_params[5], &calcdata.input_params[6]);	//read in the carbon
		fscanf(input_file, " %*[^0123456789] %lf %lf %lf\n", &calcdata.input_params[7], &calcdata.input_params[8], &calcdata.input_params[9]);	//read in the oxygen
		fprintf(swd, "\the\t%lf\t%lf\t%lf\n",calcdata.input_params[1], calcdata.input_params[2], calcdata.input_params[3]);
		fprintf(swd, "\tc \t%lf\t%lf\t%lf\n", calcdata.input_params[4], calcdata.input_params[5], calcdata.input_params[6]);
		fprintf(swd, "\to \t%lf\t%lf\t%lf\n", calcdata.input_params[7], calcdata.input_params[8], calcdata.input_params[9]);
		fprintf(swd, "\n");
		fclose(swd);
		printf("%s\n", calcdata.str_input_param.c_str());
		fscanf(input_file, "\n");
	
		//now read in desired physical properties of star
		//MUST be specified as mass, then effective temperature
		double temp;
		//read in first physical parameter -- name as string, value as double
		fscanf(input_file, "Params: %s %lf ", input_buffer, &temp);
		printf("%s %lf\n", input_buffer, temp);
		instring=std::string(input_buffer);
		//save value in appropriate slot use bit-masking to keep track of variables
		if(!instring.compare("mass")){ calcdata.mass = temp; calcdata.params|=units::pmass;}
		else {
			printf("ERROR IN INPUT: first parameter to SimpleWD must be mass\n");
			printf("This error is fatal.  Quitting.\n");
			return 1;
		}
		//read in second physical parameters -- name as string, value as double
		fscanf(input_file, "%s %lf\n", input_buffer, &temp);
		instring=std::string(input_buffer);
		//save valeu in appropriate slot, again use bit-masking so that binary value of params indicates which were used
		if(!instring.compare("teff")){ calcdata.teff = temp; calcdata.params|=units::pteff;}
		else {
			printf("ERROR IN INPUT: second parameter to SimpleWD must be teff\n");
			printf("This error is recoverable.  Using teff = 12 000 K\n");
			calcdata.teff = 12000.0;
			calcdata.params|=units::pteff;
		}
		//double check we have both mass and effective temperature
		char tempparam = units::pmass|units::pteff;
		if(calcdata.params != tempparam){
			printf("ERROR IN INPUT: a SimpleWD needs both mass and teff\n");
			printf("This error is fatal.  Quitting.\n");
			return 1;	
		}
		printf("Making SimpleWD model with:\n");
		printf("\tM=%le\n", calcdata.mass);
		printf("\tTeff=%le\n", calcdata.teff);		
	}
		
	//read the units to be used in calculation
	fscanf(input_file, "Units: %s\n", input_buffer);
	instring = std::string(input_buffer);
	if(!instring.compare("astro"))    calcdata.units = units::astro;
	else if(!instring.compare("geo")) calcdata.units = units::geo;
	else if(!instring.compare("SI"))  calcdata.units = units::SI;
	else if(!instring.compare("CGS")) calcdata.units = units::CGS;
	else{
		printf("ERROR IN INPUT: units incorrectly written...\n");
		printf("units: %s\n", input_buffer);
		printf("units: %s\n", instring.c_str());
		printf("This error is fatal.  Quitting.\n");
		return 1;
	}

	//read in the type of mode calculation to be performed
	fscanf(input_file, "Frequencies: %s", input_buffer);
	instring = std::string(input_buffer);
	if(!instring.compare("radial"))           calcdata.modetype = modetype::radial;
	else if(!instring.compare("nonradial"))   calcdata.modetype = modetype::nonradial;
	else if(!instring.compare("cowling"))     calcdata.modetype = modetype::cowling;
	else{
		printf("ERROR IN INPUT: THRAIN can only compute nonradial modes.\n");
		printf("This error is recoverable.  Changing mode type.");
		calcdata.modetype = modetype::nonradial;
	}
	//read in the adiabatic index for the mode calculation
	fscanf(input_file, "%s\n", input_buffer);
	int numerator=0,denominator=0;
	sscanf(input_buffer, "%d/%d", &numerator,&denominator);
	if(denominator!=0) calcdata.adiabatic_index = double(numerator)/double(denominator);
	else sscanf(input_buffer, "%lf", &calcdata.adiabatic_index);
	
	//Now handle reading the list of frequencies
	int startoflist = ftell(input_file);	//save location at start of list
	calcdata.mode_num = 0;					//count the number of modes to include
	while(!feof(input_file)){
		fscanf(input_file, "%s\n", input_buffer);
		calcdata.mode_num++;
	}
	printf("Number of frequencies: %d\n", calcdata.mode_num);
	fseek(input_file, startoflist, SEEK_SET); //return to start of list
	calcdata.l = new int[calcdata.mode_num];  //create arrays for mode numbers, L,K
	calcdata.k = new int[calcdata.mode_num];
	//now read in the L,K as specified 
	// calcdata.l is a pointer so addresses are passed as calcdata.l+j
	for(int j=0; j<calcdata.mode_num; j++){
		fscanf(input_file, "%d,%d\n", calcdata.l+j, calcdata.k+j);
	}

	//we are now finished reading in the input
	fclose(input_file);
	return 0;
}

//this function will print back the input file, so the same calculation can be run again
int echo_input(CalculationInputData &calcdata){
	printf("Copying input file...\t"); fflush(stdout);
	//open file to write output summary
	char output_file_name[128];
	sprintf(output_file_name, "./output/%s/%s_in.txt", calcdata.calcname.c_str(),calcdata.calcname.c_str());
	FILE* output_file;
	
	//try to open the output file
	if( !(output_file = fopen(output_file_name, "w")) ){
		//if an error occurs, try making the folder needed
		char command[140]; 
		sprintf(command, "mkdir -p ./output/%s", calcdata.calcname.c_str());
		system(command);
		if( !(output_file = fopen(output_file_name, "w")) ){
			printf("output file not found.\n");
			return 1;
		}
	}
	fprintf(output_file, "Name:\t%s\n", calcdata.calcname.c_str());
	fprintf(output_file, "Model:\t");
	switch(calcdata.regime){
		case regime::PN0:
			fprintf(output_file, "newtonian ");
			break;
	}
	switch(calcdata.model){
		case model::polytrope:
			fprintf(output_file, "polytrope %0.2lf ", calcdata.input_params[0]);
			break;
		case model::CHWD:
			fprintf(output_file, "CHWD %lf %d ", calcdata.input_params[0], (int)calcdata.input_params[1]);
			break;
		case model::MESA:
			fprintf(output_file, "MESA %s ", calcdata.str_input_param.c_str());
			break;
		case model::SWD:
			fprintf(output_file, "SWD ");
	}
	fprintf(output_file, "%d\n", calcdata.Ngrid);
	
	switch(calcdata.model){
		case model::polytrope:
			fprintf(output_file, "Params:");
			if(calcdata.params&0b0001) fprintf(output_file, "\tmass %lg", calcdata.mass);
			if(calcdata.params&0b0010) fprintf(output_file, "\tradius %lg", calcdata.radius);
			if(calcdata.params&0b0100) fprintf(output_file, "\tzsurf %lg", calcdata.zsurf);
			if(calcdata.params&0b1000) fprintf(output_file, "\tlogg %lg", calcdata.logg);
			break;
		case model::SWD:
			fprintf(output_file, "EOS:\n");
			char line1[258], line2[258];
			sscanf(calcdata.str_input_param.c_str(), "%[^\n]\n%[^\n]", line1, line2); 
			//fprintf(output_file, "%s\n", calcdata.str_input_param.c_str());
			fprintf(output_file, "\tcore\t%s\n\tatm\t%s\n", line1, line2);
			fprintf(output_file, "chemical parameters:\n");
			fprintf(output_file, "\the\t%lg\t%lg\t%lg\n",calcdata.input_params[1], calcdata.input_params[2], calcdata.input_params[3]);
			fprintf(output_file, "\tc\t%lg\t%lg\t%lg\n", calcdata.input_params[4], calcdata.input_params[5], calcdata.input_params[6]);
			fprintf(output_file, "\to\t%lg\t%lg\t%lg\n", calcdata.input_params[7], calcdata.input_params[8], calcdata.input_params[9]);
			fprintf(output_file, "Params:");
			fprintf(output_file, "\tmass %lg", calcdata.mass);
			fprintf(output_file, "\tteff %lg", calcdata.teff);
			break;
		default: break;
	}
	
	fprintf(output_file, "\nUnits:\t");
	switch(calcdata.units){
		case units::astro:
			fprintf(output_file, "astro\n");
			break;
		case units::geo:
			fprintf(output_file, "geo\n");
			break;
		case units::CGS:
			fprintf(output_file, "CGS\n");
			break;
		case units::SI:
			fprintf(output_file, "SI\n");
			break;
	}
	fprintf(output_file, "\n");
	fprintf(output_file, "Frequencies:\t");
	switch(calcdata.modetype){
		case modetype::radial:
			fprintf(output_file, "radial\t");
			break;
		case modetype::nonradial:
			fprintf(output_file, "nonradial\t");
			break;
		case modetype::cowling:
			fprintf(output_file, "cowling\t");
			break;	
	}
	fprintf(output_file, "%1.3lf\n", calcdata.adiabatic_index);
	for(int j=0; j<calcdata.mode_num; j++){
		fprintf(output_file, "%d,%d\n", calcdata.l[j], calcdata.k[j]);
	}
	
	printf("done\n");
	fflush(output_file);
	fclose(output_file);
	return 0;
}

int setup_output(CalculationInputData &data_in, CalculationOutputData &data_out){
	printf("Preparing calculation data...\n"); fflush(stdout);
	//read in the basic properties for the calculation
	data_out.calcname = data_in.calcname;
	data_out.regime = data_in.regime;
	data_out.model = data_in.model;
	data_out.modetype = data_in.modetype;
	data_out.units = data_in.units;
	data_out.num_input_params = data_in.num_input_params;
	data_out.input_params = new double[data_out.num_input_params];
	for(int a=0; a<data_out.num_input_params; a++){
		data_out.input_params[a] = data_in.input_params[a];
	}
	data_out.str_input_param = data_in.str_input_param;
	//data_out.index = data_in.index;
	data_out.Ngrid = data_in.Ngrid;
	
	//set observable parameters and format the units
	data_out.mass = data_in.mass;
	data_out.radius = data_in.radius;
	data_out.zsurf = data_in.zsurf;
	data_out.logg = data_in.logg;
	data_out.teff = data_in.teff;
	data_out.params = data_in.params;
	//formatting units may need to be re-performed after calculation, depending on star
	format_units(data_out);
	
	//now prepare the modes
	data_out.mode_num = data_in.mode_num;
	data_out.adiabatic_index = data_in.adiabatic_index;
	data_out.mode_done = 0;
	data_out.mode_writ = 0;
	data_out.mode.resize(data_out.mode_num);
	data_out.l = new int[data_out.mode_num];
	data_out.k = new int[data_out.mode_num];
	data_out.w = new double[data_out.mode_num];
	data_out.f = new double[data_out.mode_num];
	data_out.period = new double[data_out.mode_num];
	data_out.mode_SSR = new double[data_out.mode_num];
	for(int j=0; j<data_out.mode_num; j++){
		data_out.l[j] = data_in.l[j];
		data_out.k[j] = data_in.k[j];
	}
	
	//setup error columns
	data_out.i_err = 0;
	//if the star is a simple model, use RMSR to estimate mode error
	data_out.error[error::isRMSR] = ((data_out.model==model::polytrope) | (data_out.model==model::CHWD));
	//if the star is a realistic model, use overlap c_0 to estimate mode error
	data_out.error[error::isC0   ] = ((data_out.model==model::MESA) | (data_out.model==model::SWD));
	//if it is a polytrope with n=0, use the Pekeris formula to compare
	data_out.error[error::isIsopycnic] = ((data_out.model==model::polytrope) & (data_out.input_params[0]==0.0));
	//if it is a Newtonian polytrope with Gamma=5/3 and n=1.5,3,4, then compare to JCD-DJM
	data_out.error[error::isJCD] = (data_out.model==model::polytrope) &
		(data_out.regime==regime::PN0) &
		(data_out.input_params[0]==1.5 | data_out.input_params[0]==3.0 | data_out.input_params[0]==4.0) &
		(fabs(data_out.adiabatic_index - 5./3.)<1.e-5);

	//count the number of pertinent errors
	for(int e=0; e<error::numerror; e++)
		if(data_out.error[e]) data_out.i_err++;
	//create the error columns
	data_out.err = new double*[data_out.i_err];
	for(int e=0; e<data_out.i_err; e++) data_out.err[e] = new double[data_out.mode_num];
	
	printf("done\n");
	return 0;
}

int write_output(CalculationOutputData &calcdata){
	int stat=0;
	if(write_stellar_output(calcdata)) stat++;
	if(write_mode_output(calcdata)) stat++;
	return stat;
}


char dwarf[12][27] = {
//use raw stirng literal input format R"()" to make ASCII art align
	R"(#|*     * *     *|#)",
	R"(#| *    ___    * |#)",
	R"(#|  (\/.....\/)  |#)",
	R"(#|  \/.......\/  |#)", 
	R"(#| * |.-----.| * |#)", 
	R"(#|  (  U | X  )  |#)",
	R"(#|* (    v    ) *|#)",
	R"(#|_ /\_//-\\_/\ _|#)",
	R"(#|%|  /     \  |%|#)", 
	R"(#|%|     |     |%|#)",
	R"(#|%/___\|O|/___\%|#)",
	R"(#|%%%%%%%V%%%%%%%|#)"
};

int write_stellar_output(CalculationOutputData& calcdata){
	int d=0;
	printf("Writing stellar data to file...\t");fflush(stdout);
	//open file to write output summary
	char output_file_name[128];
	sprintf(output_file_name, "./output/%s/%s.txt", calcdata.calcname.c_str(),calcdata.calcname.c_str());
	FILE* output_file;
	//try to open the output file
	if( !(output_file = fopen(output_file_name, "w")) ){
		//if an error occurs, try making the folder needed
		printf("creating file..."); fflush(stdout);
		char command[140]; sprintf(command, "mkdir -p ./output/%s", calcdata.calcname.c_str());
		system(command);
		if( !(output_file = fopen(output_file_name, "w")) ){
			printf("output file not found.\n");
			return 1;
		}
	}
	
	//print the title bar
	int WIDTH = 86;
	//first bar
	fprintf(output_file, "#");
	for(int j=1; j<WIDTH; j++) fprintf(output_file, "-");
	//code title
	std::string splash = "|    THRAIN     |  The Mighty White Dwarf Code";
	fprintf(output_file, "\n#%s \n", splash.c_str());
	//second bar
	fprintf(output_file, "#");
	for(int j=1; j<WIDTH; j++) fprintf(output_file, "-");
	fprintf(output_file, "\n");
	
	char s[256],m[256];
	switch(calcdata.model){
		case model::polytrope:
		sprintf(s,"Polytropic star");
		break;
		case model::CHWD:
		sprintf(s,"Chandrasekhar WD");
		break;
		case model::MESA:
		sprintf(s,"MESA model");
		break;
		case model::SWD:
		sprintf(s,"Simple WD model");
		break;
	}
	switch(calcdata.modetype){
		case modetype::radial:
			sprintf(m,"radial");
			break;
		case modetype::nonradial:
			sprintf(m,"nonradial");
			break;
		case modetype::cowling:
			sprintf(m,"Cowling");
			break;
	}	
	fprintf(output_file, "%s  %s with %s pulsations\n", dwarf[d++], s, m);

	//print out a message showing which observable parameters were passed, and with which units
	char unitM[100], unitL[10], unitT[10], unitG[10], unitZ[10];
	switch(calcdata.units){
		case units::astro:
			sprintf(unitM, "(Msolar)");
			sprintf(unitL, "(km)    ");
			sprintf(unitG, "(cm/s^2)");
			sprintf(unitZ, "        ");
			//fprintf(output_file,"%s  Astronomical Units (Msolar, km, s)\n", dwarf[d++]);
			break;
		case units::geo:
			sprintf(unitM, "(m)     ");
			sprintf(unitL, "(m)     ");
			sprintf(unitG, "(cm/s^2)");
			sprintf(unitZ, "        ");
			//fprintf(output_file,"%s  Geometric Units (G=c=1)\n", dwarf[d++]);
			break;
		case units::SI:
			sprintf(unitM, "(kg)    ");
			sprintf(unitL, "(m)     ");
			sprintf(unitG, "(cm/s^2)");
			sprintf(unitZ, "        ");
			//fprintf(output_file,"%s  SI Units (kg, m, s)\n", dwarf[d++]);
			break;
		case units::CGS:
			sprintf(unitM, "(g)     ");
			sprintf(unitL, "(cm)    ");
			sprintf(unitG, "(cm/s^2)");
			sprintf(unitZ, "        ");
			//fprintf(output_file,"%s  CGS Units (g, cm, s)\n", dwarf[d++]);
			break;
	}
	fprintf(output_file,"%s      Mass   %s = %1.5lg %s", dwarf[d++], unitM, calcdata.mass,  (calcdata.params&units::pmass?"(specified)\n":"(derived)\n"));
	fprintf(output_file,"%s      Radius %s = %1.5lg %s", dwarf[d++], unitL, calcdata.radius,(calcdata.params&units::pradius?"(specified)\n":"(derived)\n"));
	if(calcdata.teff!=0.0)
	fprintf(output_file,"%s      Teff (K)%s= %lg %s",    dwarf[d++], unitZ, calcdata.teff,  (calcdata.params&units::pteff?"(specified)\n":"(derived)\n"));
	else
	fprintf(output_file,"%s      Teff      = N/A \n",    dwarf[d++]);	
	fprintf(output_file,"%s      log g%s   = %1.5lg %s", dwarf[d++], unitG, calcdata.logg,  (calcdata.params&units::plogg?"(specified)\n":"(derived)\n"));
	fprintf(output_file,"%s      Zsurf  %s = %1.5le %s", dwarf[d++], unitZ, calcdata.zsurf, (calcdata.params&units::pzsurf?"(specified)\n":"(derived)\n"));
	fprintf(output_file,"%s  \n", dwarf[d++]);
	

	//printf the background model data
	switch(calcdata.model){
		case model::polytrope:
			fprintf(output_file, "%s  index n = %.2lf\n", dwarf[d++], calcdata.input_params[0]);
			fprintf(output_file, "%s  \n", dwarf[d++]);
			break;
		case model::CHWD:
			fprintf(output_file, "%s  y0 = %.2lf", dwarf[d++], calcdata.input_params[0]);
			fprintf(output_file, "%s  \n", dwarf[d++]);
			break;
		case model::MESA:
			fprintf(output_file, "%s  model file %s.dat", dwarf[d++], calcdata.str_input_param.c_str());
			fprintf(output_file, "%s  \n", dwarf[d++]);
			break;
		case model::SWD:
			fprintf(output_file, "%s  Layer masses: \n", dwarf[d++]);
			fprintf(output_file, "%s    ", dwarf[d++]);
			fprintf(output_file, "H: %1.3le He: %1.3le  C: %1.3le  O: %1.3le\n", 
				((SimpleWD*) calcdata.star)->Xmass.H1,
				((SimpleWD*) calcdata.star)->Xmass.He4,
				((SimpleWD*) calcdata.star)->Xmass.C12,
				((SimpleWD*) calcdata.star)->Xmass.O16);
			break;
	}
	fprintf(output_file, "%s  \n", dwarf[d++]);
	fprintf(output_file, "%s  Number of grid points : %d\n", dwarf[d++], calcdata.Ngrid);
	fprintf(output_file, "%s  Fractional RMS error  : %1.3le\n", dwarf[d++], calcdata.star_SSR);
	//bottom bar
	fprintf(output_file, "#");
	for(int j=1; j<WIDTH; j++) fprintf(output_file, "-");
	fprintf(output_file, "\n");
	fclose(output_file);
	
	char outname[128];
	sprintf(outname, "./output/%s",calcdata.calcname.c_str());
	calcdata.star->writeStar(outname);
	
	printf("done\n");
	return 0;
}

int write_mode_output(CalculationOutputData& calcdata){
	printf("Writing mode data to file...\t"); fflush(stdout);
	//open file to write output summary
	char output_file_name[128];
	sprintf(output_file_name, "./output/%s/%s.txt", calcdata.calcname.c_str(),calcdata.calcname.c_str());
	FILE* output_file;
	//try to open the output file
	if( !(output_file = fopen(output_file_name, "a")) ){
		printf("the file doesn't exist\n");
		return 1;
	}
		
	int WIDTH = 86;
	if(calcdata.mode_writ==0){
		fprintf(output_file, "#  Stellar Pulsation Results  \n");
		//Printf the adiabatic index used in calculation -- 5/3, 4/3 special cases
		fprintf(output_file, "#  Adiabatic exponent (Gamma1) ");
		double a = calcdata.adiabatic_index*3.0;
		if(fabs(a-5.0)<1e-10) fprintf(output_file, "= 5/3\n");
		else if (fabs(a-4.0)<1e-10) fprintf(output_file, "= 4/3\n");
		else if (a==0.0) fprintf(output_file, "matched to stellar profile\n");
		else fprintf(output_file, "= %lf\n", calcdata.adiabatic_index);
		fprintf(output_file, "#");
		for(int j=1; j<WIDTH; j++) fprintf(output_file, "-");
		fprintf(output_file, "\n");
		
		//begin printing titles for columns
		//line 1
		fprintf(output_file, "#     \tmode \tfreq               \tfreq            \tperiod             ");
		std::vector<std::string> topline{
			"\tfractional  ",
			"\t            ",
			"\trel.error w/",
			"\tabs.diff. w/"
		};
		for(int e=0; e<error::numerror; e++)
			if(calcdata.error[e]) fprintf(output_file, "%s", topline[e].c_str());                                                                 
		fprintf(output_file, "\n");                                                                                    
		//line 2                                                                                                        
		fprintf(output_file, "# L,N \ttype \tsq(GM/R^3)         \t(Hz)            \t(s)               ");
		std::vector<std::string> botline{
			"\tRMSR     ",
			"\tc0 or c1 ",
			"\tPekeris formula",
			"\tJCD_DJM 1994 (uHz)"
		};
		for(int e=0; e<error::numerror; e++)
			if(calcdata.error[e]) fprintf(output_file, "%s", botline[e].c_str());
		fprintf(output_file, "\n#");
		for(int j=1; j<WIDTH; j++) fprintf(output_file, "-");
		fprintf(output_file, "\n");
	}
	
	int start = calcdata.mode_writ;
	for(int j=start; j<calcdata.mode_done; j++){
		//print the mode numbers L,K
		fprintf(output_file, " %d,%d \t", calcdata.l[j], calcdata.k[j]);
		//print a mode-type label (p,f,g)
		     if (calcdata.k[j] <0.0) fprintf(output_file, "g   \t");
		else if (calcdata.k[j] >0.0) fprintf(output_file, "p   \t");
		else if (calcdata.k[j]==0.0) fprintf(output_file, "f   \t");
		if(calcdata.w[j]==0.0){
			fprintf(output_file, "unable to find mode\n");
			calcdata.mode_writ++;
			continue;
		}
		fprintf(output_file, "%0.12le \t%3.12le \t%0.12le", calcdata.w[j], calcdata.f[j],calcdata.period[j]);
		//fprintf(output_file, "\t%1.2le", calcdata.mode_SSR[j]);
		for(int e=0; e<calcdata.i_err; e++){
			if(!std::isnan(calcdata.err[e][j])) fprintf(output_file, "\t%1.2le", calcdata.err[e][j]);
			else fprintf(output_file, "\tN/A");
		}
		fprintf(output_file, "\n");
		char outname[128];
		sprintf(outname, "./output/%s",calcdata.calcname.c_str());
		(calcdata.mode[j])->writeMode(outname);
		fflush(output_file);
		calcdata.mode_writ++;
	}

	fclose(output_file);
	printf("done\n");
	return 0;
}

void print_splash(FILE* output_file, char* title, int WIDTH){
	//first bar
	fprintf(output_file, "#");
	for(int j=1; j<WIDTH; j++) fprintf(output_file, "-");
	//code title
	std::string splash = "|    THRAIN     |  The Mighty White Dwarf Code";
	fprintf(output_file, "\n#%s \n", splash.c_str());
	//second bar
	fprintf(output_file, "#");
	for(int j=1; j<WIDTH; j++) fprintf(output_file, "-");
	fprintf(output_file, "\n");
}


//NOTE: this function has a problem and will always cause a seg fault
int write_tidal_overlap(CalculationOutputData& calcdata){
	printf("Writing tidal overlap coefficients...\n");fflush(stdout);
	//open file to write output summary
	char output_file_name[128];
	sprintf(output_file_name, "./output/%s/tidal_overlap.txt", calcdata.calcname.c_str());
	FILE* output_file;
	//try to open the output file
	if( !(output_file = fopen(output_file_name, "w")) ){
		//if an error occurs, try making the folder needed
		printf("creating file...\n"); fflush(stdout);
		char command[140]; sprintf(command, "mkdir -p ./output/%s", calcdata.calcname.c_str());
		system(command);
		if( !(output_file = fopen(output_file_name, "w")) ){
			printf("output file not found.\n");
			return 1;
		}
	}
	
	//print the cool splash
	int WIDTH = 80;
	char splashy[1100];
	sprintf(splashy, "tidal overlap coefficients for %s", calcdata.calcname.c_str());
	print_splash(output_file, splashy, WIDTH);
	fflush(output_file);
		
	//a useful error measurement is to calculate c0 -- see Fuller & Lai 2011 A
	//for this, we need the k=0 mode for each value of l
	//produce a list of the different L asked for
	int lastl=0, minl=100, maxl=-1;
	//find min, max
	int l_current = calcdata.l[0];
	if(l_current < minl) minl=l_current;
	if(l_current > maxl) maxl=l_current;
	int num = calcdata.mode_done;
	int l_list[num];
	printf("\tpreparing mode list...\t");
	l_list[0] = calcdata.l[0];
	for(int j=1; j<num; j++){
		int l_current = calcdata.l[j];
		if(l_current == l_list[lastl]) continue;
		// do not calculate overlap for dipole or radial modes
		if(l_current < 2) {l_list[lastl]=l_current; continue;}
		else{
			//check if we already counted this l value
			bool already = false;
			for(int i=0; (i<=lastl) & !already; i++){
				already |= (l_current == l_list[i]);
			}
			if(already) continue;
			else {
				lastl++;
				l_list[lastl] = l_current;
			}
		}
		//find min, max
		if(l_current < minl) minl=l_current;
		if(l_current > maxl) maxl=l_current;
	}
	//extra increment to align with zero-indexing
	lastl++;
	for(int j=lastl; j<num; j++) l_list[j]=-1;
	printf("done\n");
	
	//an array to allow indexing of fmodes by l
	int jforl[maxl-minl+1];
	for(int l=minl; l<=maxl; l++){
		for(int j=0; j<lastl; j++){
			if(l==l_list[j]) jforl[l] = j;
		}
	}
		
	//now make the fundamental mode for each of those l
	ModeBase *fmode[lastl];
	printf("\tpreparing f-modes...\n");fflush(stdout);
	for(int j=0; j<lastl; j++){
		printf("\t\tl=%d\t", l_list[j]);
		//for dipole (l=1) modes, there is no f-mode-- skip!
		if(l_list[j]<2) {
			fmode[j] = NULL;
			printf("no tidal response at this order\n");
			continue;
		}
		//check if the f-mode has already been found earlier
		bool inlist=false;
		for(int i=0; i<num; i++){
			if(calcdata.l[i]==l_list[j] & calcdata.k[i]==0) {
				inlist=true;
				fmode[j] = calcdata.mode[i];	
			}
		}
		//if so, we don't need to calculate it, move on to next l
		if(inlist) {
			printf("already found!\n");
			continue;
		}
		//if not, then keep searching until we find it
		printf("searching ...\t");fflush(stdout);

		fmode[j] = new Mode<4>(0,l_list[j],0, calcdata.driver);	
		int kk = fmode[j]->modeOrder();
		if(kk!=0){
			//use a bisection search to find the desired mode
			double w2min=0.0, w2max=0.0, dw2=0.0, w2in=0.0, w2out=0.0;
			int kmax=0, kmin=0;
			//first we create brackets
			//bracket search when discovered mode is HIGHER than desired mode
			if(kk > 0){		
				//use the current mode as a max bracket (since it is high)
				kmax  = kk;
				w2max = fmode[j]->getOmega2();
				//use an absolute minimum as a min bracket
				// (note: this puts limits on allowed gmodes at -999999999)
				w2min = 0.0;
				kmin = -1000000000;
			}
			//bracket search when discovered mode is LOWER than desired mode
			else if(kk < 0){
				//use the current mode as a min bracket
				kmin  = kk;
				w2min = fmode[j]->getOmega2();
				dw2 = w2min;
				//if the max bracket was not found in list, search for it
				//start at current mode and increase
				kmax  = kk;
				w2max = w2min+dw2;
				//increase and search until we find a max bracket
				while(kmax < 0){
					delete fmode[l_list[j]];
					w2max = w2max + dw2;

					//fmode[l_list[j]] = new Mode<NV1>(w2max, l_list[j],0,calcdata.driver);
					kmax = fmode[j]->modeOrder();
					//if we found it, quit
					if(kmax == 0){
						kk=kmax;
						break;
					}
				}
			}
			if(kmax==0) break;
			//swap brackets if backward
			if(kmin >= kmax){
				int tk = kmax;
				kmax = kmin;
				kmin = tk;
				double tw = w2max;
				w2max = w2min;
				w2min = tw;
			}
		
			//now we have brackets -- these SHOULD put bounds in frequency
			//for w2 in (w2min, w2max), will produce k in (kmin, kmax)
			double prevmin=w2min, prevmax=w2max; //value of previous frequency
			int stop=0; //integer to limit number of iterations
			while(kk != 0){
				w2in = 0.5*(w2min+w2max); //bisect the brackets
				//create a trial mode
				delete fmode[j];
				fmode[j] = new Mode<4>(w2in, l_list[j],0,calcdata.driver);
				kk = fmode[j]->modeOrder();
				w2out = fmode[j]->getOmega2();
			
				//if we found it, then great.  move on to next
				if(kk == 0) {
					break;
				}
				//if we didn't find it, see if either bracket can be moved
				if(     kk > 0 & kk <= kmax & w2max >w2out & w2out>0.0){
					kmax = kk;
					w2max = w2out;
				}
				else if(kk < 0 & kk >= kmin & w2min<w2out & w2out>0.0){
					kmin = kk;
					w2min = w2out;
				}
											
				//if we are just unable to find the mode, say so
				if(fabs(w2max-w2min) < 1e-10){
					printf("too close\t%le\n", fabs((w2max-w2min)/w2max) );
					break;
				}
				//update past values
				prevmin = w2min, prevmax=w2max;
			}
		}
	}
	printf("\tdone\n");
	
	printf("\tcalculating overlap and c0...\t");
	fprintf(output_file, "#l,k \tmodeid\tomega^2 (GM/r^3)  \tdimensionless overlap \tc0\n");
	for(int j=0; j<WIDTH; j++) fprintf(output_file, "#");
	fprintf(output_file, "\n");fflush(output_file);
	int ll=l_list[0];
	for(int j=0; j<calcdata.mode_done; j++){
		if(calcdata.l[j]<2) continue;
		//print the mode numbers L,K
		fprintf(output_file, " %d,%d \t", calcdata.l[j], calcdata.k[j]);
		//print a modetype label (p,f,g)
		     if (calcdata.k[j]<0.0)  fprintf(output_file, "g    \t");
		else if (calcdata.k[j]>0.0)  fprintf(output_file, "p    \t");
		else if (calcdata.k[j]==0.0) fprintf(output_file, "f    \t");
		if(calcdata.w[j]==0.0){
			fprintf(output_file, "unable to find mode\n");
			continue;
		}
		fprintf(output_file, "%0.12le \t", sqrt(calcdata.mode[j]->getOmega2())); fflush(output_file);
		fprintf(output_file, "%3.16le \t", fabs(calcdata.mode[j]->tidal_overlap())); fflush(output_file);
		fprintf(output_file, "%0.12le \n", 
			fabs(calcdata.driver->innerproduct(calcdata.mode[j],fmode[jforl[calcdata.l[j]]]))
		);
		fflush(output_file);
	}
	fclose(output_file);
	for(int j=0; j<lastl; j++){
		delete fmode[j];
	}
	printf("\tdone\n");
	printf("done!\n");
	return 0;
}
