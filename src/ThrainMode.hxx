//**************************************************************************************
// THE MODE FINDER ALGORITHM	
// This function is the true workhorse of THRAIN.
// This will organize the specified mode numbers, then begin looking through the list
// If it does not first find the desired mode, it will use a bisection search to find it
// If an accidentally discovered mode is in a list, it will save it
// Mode output is printed at the end of each run of constant L
//**************************************************************************************

#ifndef THRAINMODEHXX
#define THRAINMODEHXX

#include "ThrainMain.h"

namespace rootfind {
	double pseudo_unif();
}

namespace mode {

template <class MODE>
void save_mode (
	MODE* modetry,
	std::set<int>& kfilled,
	std::map<int,double>& w2filled,
	std::map<int, MODE*>& modefilled
){
	static_assert((std::is_base_of<ModeBase,MODE>::value), "first class must be mode");
	int K = modetry->modeOrder();
	double w2 = modetry->getOmega2();
	// add to the lists
	if(!kfilled.count(K) && w2>0.0){ // zero frequencies represent an error
		kfilled.insert(K);
		w2filled[K] = w2;
		modefilled[K] = modetry;
		modetry = nullptr;
	}
	delete modetry;
}

template <class MODE, class MODEDRIVER>
int mode_finder(Calculation::OutputData &data){
	//ensure that the passed classes have correct parents and are compatible
	static_assert((std::is_base_of<ModeDriver,MODEDRIVER>::value), "second class must be modedriver");
	static_assert((std::is_base_of<ModeBase,MODE>::value), "first class must be mode");
	static_assert((MODE::num_var==MODEDRIVER::num_var), "incompatible mode, mode class");

	printf("num_var = %d\n", MODE::num_var);
	
	//information on the nodes that need to be found
	int num = data.mode_num;
	// good ol' two pi
	double twopi = 6.283185307179586;
	
// STEP 1: clean up the input
// produce a list of the different L,K asked for, each represented once
	std::set<int> l_list;
	std::map<int, std::set<int>> kl_lists;
	for(int j=0; j<num; j++){
		l_list.insert(data.l[j]);
		kl_lists[data.l[j]].insert(data.k[j]);
	}
	// always add an f-mode, for error estimation
	for(auto lt=l_list.begin(); lt!=l_list.end(); lt++){
		kl_lists[*lt].insert(0);
	}
	// remove the dipole f-mode -- use p1 mode instead
	if(l_list.count(1)){
		kl_lists[1].erase(0);
		kl_lists[1].insert(1);
	}
	// remove the list elements
	data.l.clear();
	data.k.clear();
	data.l.reserve(data.mode_num);
	data.k.reserve(data.mode_num);

	//now for each L, find the indicated K
	int nextmode = 0;
	for(auto lt=l_list.begin(); lt!=l_list.end(); lt++){
		const int ltarget = *lt;
		printf("L=%d\n", ltarget);

		std::set<int> kl = kl_lists[ltarget];
		//these will be filled as we discover modes
		// each is mapped to by its k-value 
		std::set<int> kfilled;
		std::map<int,double> w2filled;
		std::map<int, MODE*> modefilled;

		int ktry;
		double w2try;
		MODE* modetry = nullptr;
		
// STEP 2:  perform first run with rough guesses 
// STEP 2a: fill in all easy modes from simple guesses based on K
		printf("\tinitial calculation of modes...\t"); fflush(stdout);
		for(auto kt=kl.begin(); kt!=kl.end(); kt++){
			modetry = new MODE(*kt, ltarget, 0, data.driver);
			mode::save_mode<MODE>(modetry, kfilled, w2filled, modefilled);
		}
		std::size_t found=0, total=0;
		for(auto kt=kl.begin(); kt!=kl.end(); kt++){
			total++;
			found += kfilled.count(*kt);
		}
		printf("found %lu/%lu!\n", found, total);

		// save values for the min,max values of K in list
		int klo=*(kl.begin()), khi=*(kl.rbegin());


		//****************************************************************************
		// This will run through each k sublist, attempting new modes
		//   The calculation does not always return the desired mode
		//   However, the returned mode is often in the list of modes
		//   This adjusts the search for correct k, while also saving those found by
		//     accident to reduce the amount of recalculation
		//   Search uses a bracketed bisection search on omega2.
		//****************************************************************************

		for(auto kt=kl.begin(); kt!=kl.end(); kt++){	
			int ktarget = *kt;
			printf("%d\tK=%d\t", ltarget, ktarget); fflush(stdout);

// STEP 2b:  check if we have already been filled from earlier
			if(kfilled.count(ktarget)){
				//if we have then continue to next
				printf("k=%d\talready found!\n", ktarget);
				continue;
			}
		
// STEP 3: set brackets
// STEP 3a: search within discovered modes to see if brackets exist	
			printf("finding brackets.\n");
			ktry = klo-1;
			double w2min=0.0, w2max=0.0, dw2=0.0, w2in=0.0;
			int kmax=klo-1, kmin=khi+1;
			mode::get_min_from_set(kfilled, w2filled, ktarget, kmin, w2min);
			mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
			printf("\t\t(%d,%d) first fit\n", kmin, kmax);

// STEP 3b: if one or both brackets do not exist, look for brackets
			bool nomax = (kmax==klo-1), nomin = (kmin==khi+1); 
			if( nomax || nomin ) {
				// create a trial mode and add it to the list
				modetry = new MODE(ktarget, ltarget, 0, data.driver);
				ktry = modetry->modeOrder();
				mode::save_mode<MODE>(modetry, kfilled, w2filled, modefilled);
				//if this is the one we want, continue to next mode
				if(kfilled.count(ktarget)){
					printf("\tFOUND k=%d, w2=%lf\n", ktry,w2try);
					continue;
				}
// STEP 3c: if no min bracket, try to use current, or else use absolute min
				if(nomin){
					// see if we can use current
					if(ktry < ktarget) {
						mode::get_min_from_set(kfilled, w2filled, ktarget, kmin, w2min);
					}
					// else use absolute min
					else {
						kmin = -1000000000; // THRAIN should never be used for this order g-mode
						w2min = 0.0;
					}
				}

// STEP 3d: if no max bracket, try to use current, or else quest for max
				if(nomax){
					// see if we can use current
					if(ktry > ktarget){
						mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
					}
					// else quest for it -- keep searching iteratively higher
					else {
						//start at current highest mode and increase
						kmax  = *(kfilled.rbegin());
						w2max = w2filled[kmax];
						MODE* modeMaxQuest = nullptr; // preserve modetry while on quest
						double incr = 2.0;
						while(kmax < ktarget){
							// delete modeMaxQuest;
							// increase the frequency and look
							w2max = incr*w2max;
							// create a mode with this larger frequency
							modeMaxQuest = new MODE(w2max, ltarget, 0,data.driver);
							kmax = modeMaxQuest->modeOrder();
							w2max = modeMaxQuest->getOmega2();
							printf("\t\t(%d,%d) in (%f,%f)\n",kmin, kmax, w2min, w2max);
							if( std::isnan(w2max) ) return 1;
							// add this mode to discovered lists
							mode::save_mode<MODE>(modeMaxQuest, kfilled, w2filled, modefilled);
							//if we found it, quit
							if(kmax == ktarget){
								break;
							}
							// if we exceeded bounds, stop		
							if(w2max > 1e6) {
								printf("TOO LARGE\n");
								incr = 0.5;
							}
							if(w2max < 1e-6){
								printf("No Dice\n");
								break;
							}
							mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
						}
						mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
					}
				}
			}
		
// STEP 3c: if we could not find a max bracket, go on to next mode
			if(kmax < ktarget){
				printf("\t\tunable to find max bracket\n");
				continue;
			}
			printf("\t\tbracketed (%d,%d) in (%f,%f)\n", kmin, kmax, w2min, w2max);

// STEP 4:  now we have brackets -- these SHOULD put bounds in frequency
// for w2 in (w2min, w2max), will produce k in (kmin, kmax)
// reiterate bisection search until desired mode is found
			double prevmin=w2min, prevmax=w2max; //value of previous frequency
			int stop=0; //integer to limit number of iterations
			while(!kfilled.count(ktarget)){
//STEP 4a: bisect the brackets 
				w2in = 0.5*(w2min+w2max); //bisect the brackets

//STEP 4b: create a trial mode and add to list
				modetry = new MODE(w2in, ltarget,0,data.driver);
				ktry = modetry->modeOrder();
				w2try = modetry->getOmega2();
				mode::save_mode<MODE>(modetry, kfilled, w2filled, modefilled);
				printf("%d\t\t(%d,%d) in (%f,%f)\t %d\t %d %lf-->%lf\n",
						ltarget,kmin, kmax, w2min, w2max, *kt, ktry, w2in, w2try);

//STEP 4c: compare the trial mode to desired mode
//if we found it, move on to next
				if(ktry==ktarget) {
					printf("%d\tK=%d\tk=%d FOUND\n", ltarget, ktarget, ktry);
					break;
				}

//STEP 4d: if we didn't find it, try to move brackets
				mode::get_min_from_set(kfilled, w2filled, ktarget, kmin, w2min);
				mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
				//accounts for fact multiple w2in lead to same k
				if(ktry > ktarget && w2in < w2max && w2try>0.0){
					w2max = w2in;
				}
				else if(ktry < ktarget && w2in > w2min && w2try>0.0){
					w2min = w2in;
				}

//STEP 4e: if the sought mode is bracketed between two known frequencies
// then we can use special Mode constructor
				if((kmin == ktarget-1) && (kmax == ktarget+1)){
					modetry = new MODE(w2min, w2max, ltarget,0,data.driver);
					ktry = modetry->modeOrder();
					w2try = modetry->getOmega2();
					printf("%d\t\t(%d,%d) in (%f,%f)\t %d\t %d %lf-->%lf\n",
						ltarget,kmin, kmax, w2min, w2max, *kt, ktry, w2in, w2try);
					// add this to the lists
					mode::save_mode<MODE>(modetry, kfilled, w2filled, modefilled);
				}
				
//STEP 4f: check if the brackets have moved since last iteration
				//if not, pick a random location within brackets to test
				//repeat until one of the brackets can move
				//this prevents us from getting stuck
				double w2maxT = w2max, w2minT = w2min;
				int enough=0;
				while(ktry!=ktarget && (w2min==prevmin)&&(w2max==prevmax)){
					//pick a pseudo-random place in brackets
					w2in = w2minT + rootfind::pseudo_unif()*(w2maxT-w2minT);
					modetry = new MODE(w2in, ltarget,0,data.driver);
					ktry = modetry->modeOrder();
					w2try = modetry->getOmega2();
					mode::save_mode<MODE>(modetry, kfilled, w2filled, modefilled);
					printf("\tR\t(%d,%d) in (%f,%f)\t %d\t %d %lf-->%lf\n",
							kmin, kmax, w2minT, w2maxT, ktarget, ktry, w2in, w2try);
					//check if we found mode, or if we can move brackets
					if(ktry==ktarget) break;
					mode::get_min_from_set(kfilled, w2filled, ktarget, kmin, w2min);
					mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
					
					//accounts for fact multiple w2in lead to same k
					if(w2in < w2max && ktry > ktarget && w2in>0.0){
						w2max = w2in;
					}
					else if(w2in > w2min && ktry < ktarget && w2in>0.0){
						w2min = w2in;
					}
					//sometimes zeros are inaccessible
					//move brackets and try again
					if(fabs(w2min-w2max)<1e-2*w2min) {
						w2min = w2minT;
						w2max = w2maxT;
						break;
					}
					//cancel if we have tried more than 5 random spots
					if(enough++ > 5) {
						w2min=w2minT;
						w2max=w2maxT;
						break;
					}
				}
				//update past values
				prevmin = w2min, prevmax=w2max;

				// if we found it, say so
				if(ktry == ktarget) {
					printf("%d\tK=%d\tk=%d FOUND\n", ltarget, ktarget, ktry);
					break;
				}

				//if we are just unable to find the mode, say so
				if(ktry!=ktarget && fabs(w2max-w2min) < 1e-2*w2min){
					printf("too close\t%le\n", fabs((w2max-w2min)/w2max) );
					break;
				}
				
				if(stop++ > 10) {
					break;
				}
			}
		}
		
// STEP 5: return through list to pick up any missed modes
		// this calculation follows same steps as STEP 3 but because more
		// modes have been filled, we are more likely to have good brackets
		printf("L=%d\tpicking up stragglers...\n", ltarget);
		for(auto kt= kl.begin(); kt!=kl.end(); kt++){
			const int ktarget = *kt;
			//STEP 5a:  check if we have already been filled from earlier
			if(kfilled.count(ktarget)) {
				continue;
			}
			//STEP 5b: if we have not, fill
			else {
				//STEP 5b (i): find brackets to use in bisection by scanning list of modes
				//  this pass through, we are more likely to have bracketing modes
				printf("\t\t%d\t", *kt); fflush(stdout);
				ktry = klo-1;
				double w2min=0.0, w2max=0.0, dw2=0.0, w2in=0.0;
				const int khi=*(kfilled.rbegin()), klo=*(kfilled.begin());
				int kmax=khi+1, kmin=klo-1;
				
				//if the list of previously found modes (ktry) contains a lower mode
				//then use that mode as a minimum bracket
				mode::get_min_from_set(kfilled, w2filled, ktarget, kmin, w2min);
				mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
				
				//if we don't have a lower bound from the list, use absolute minimum
				if(kmin < klo){ 
					w2min = 0.0;
					kmin = -1000000000;
				}
				//if we don't have a higher bound from calculated list, just leave
				if(kmax > khi) { printf("no brackets\n"); continue;}
				printf("(%d,%d)\t", kmin,kmax);fflush(stdout);
				
				//STEP 5b (ii): now bisect the brackets until correct mode is found
				double prevmin=w2min, prevmax=w2max; //value of previous frequency
				int stop=0; //integer to limit number of iterations
				while(!kfilled.count(ktarget)){
					w2in = 0.5*(w2min+w2max); //bisect the brackets
					//create a trial mode
					modetry = new MODE(w2in, ltarget,0,data.driver);
					ktry = modetry->modeOrder();
					w2try = modetry->getOmega2();
					// if in list
					mode::save_mode<MODE>(modetry, kfilled, w2filled, modefilled);
					//if we found it, then great.  move on to next
					if(ktry == ktarget) {
						break;
					}
					//if we didn't find it, see if either bracket can be moved
					mode::get_min_from_set(kfilled, w2filled, ktarget, kmin, w2min);
					mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
					//if the sought mode is bracketed between two known frequencies
					// then we can use special Mode constructor
					if((kmin == ktarget-1) && (kmax == ktarget+1)){
						modetry = new MODE(w2min, w2max, ltarget,0,data.driver);
						ktry = modetry->modeOrder();
						w2try = modetry->getOmega2();
						mode::save_mode<MODE>(modetry, kfilled, w2filled, modefilled);
					}
					
					//check if the brackets have moved since last time
					//if not, pick a random location within brackets to test
					//repeat until one of the brackets can move
					double w2maxT = w2max, w2minT = w2min;
					int enough=0;
					while(ktry != ktarget && (w2min==prevmin)&&(w2max==prevmax)){
						//if scanning didn't work, pick a pseudo-random place in brackets
						w2in = w2minT + rootfind::pseudo_unif()*(w2maxT-w2minT);
						modetry = new MODE(w2in, *lt,0,data.driver);
						ktry = modetry->modeOrder();
						w2try = modetry->getOmega2();
						mode::save_mode<MODE>(modetry, kfilled, w2filled, modefilled);
						//check if we found mode, or if we can move brackets
						if(ktry == ktarget) break;
						mode::get_min_from_set(kfilled, w2filled, ktarget, kmin, w2min);
						mode::get_max_from_set(kfilled, w2filled, ktarget, kmax, w2max);
						//accounts for fact multiple w2in lead to same k
						if(w2in < w2maxT && ktry == kmax & w2in>0.0){
							w2maxT = w2in;
						}
						else if(w2in > w2minT && ktry == kmin & w2in>0.0){
							w2minT = w2in;
						}
						//sometimes zeros are inaccessible
						//move brackets and try again
						if(fabs(w2minT-w2maxT)<1e-2*w2minT) {
							w2min = w2minT;
							w2max = w2maxT;
							break;
						}
						//cancel if we have tried more than 20 random spots
						if(enough++ > 5) {
							w2min=w2minT;
							w2max=w2maxT;
							break;
						}
					}						
					
					//update past values
					prevmin = w2min, prevmax=w2max;
					if(stop++ > 5) {
						break;
					}
				}
				//if we found it, say so
				if(kfilled.count(ktarget)) printf("found\n");
				else printf("not found\n");
			}
		}
		printf("\tdone\n");

//STEP 6: organize the data lists
		printf("preparing output...\t"); fflush(stdout);
		int enext = data.mode_done;
		for(auto kt=kl.begin(); kt!=kl.end(); kt++){
			const int ktarget = *kt;
			data.l.push_back(ltarget);
			data.k.push_back(ktarget);
			//STEP 5a: if a mode was found, save its data
			if(kfilled.count(ktarget)){
				double w2try = w2filled.at(ktarget);
				data.w.push_back( w2try>0 ? sqrt(w2try) : -sqrt(-w2try) );
				data.mode.push_back(modefilled.at(ktarget));
				//freq0 converts to rad/s, divide by 2pi to get frequency in Hertz
				data.f.push_back(data.freq0*data.w[nextmode]/(twopi));
				//period in seconds is 1/f
				data.period.push_back(1./data.f[nextmode]);
				data.mode_SSR.push_back(data.mode[nextmode]->SSR());
			}
			//STEP 5b: otherwise, say so and save nothing to it
			else {
				data.w.push_back(0.0);
				data.f.push_back(0.0);
				data.period.push_back(0.0);
				data.mode_SSR.push_back(1.0);	
				data.mode.push_back(nullptr);
			}
			nextmode++;
			data.mode_done++;
		}
		printf("done\n");
	
//STEP 7: calculate desired errors
		printf("calculating errors...\t"); fflush(stdout);
		int e=0;
		//STEP 7a: for simple models, use RMSR to indicate numerical error
		if(data.error[error::isRMSR]){
			for(int i=enext;i<data.mode_done; i++){
				data.err[e][i] = data.mode_SSR[i];
			}
			e++;
		}
		//STEP 7b: for realistic models, use overlap c0 to indicate numerical error
		//  this logic is not working very well yet
		if(data.error[error::isC0]){
			int testK = (*lt==1 ? 1 : 0);
			MODE *testmode;
			bool containsK=false;
			int indexK = 0;
			//check if we already have the test mode to compare against
			for(int i=enext;i<data.mode_done;i++){
				if((!containsK) & (data.k[i]==testK)) indexK=i;
				containsK |= (data.k[i]==testK);
			}
			//if so use it
			if(containsK) testmode = static_cast<MODE*>(data.mode[indexK]);
			//if not, make one
			else {
				testmode = new MODE(testK,*lt,0,data.driver);
			}
			bool correctTest = (testmode->modeOrder() == testK);
			for(int i=enext;i<data.mode_done; i++){
				if(data.w[i]==0) continue;
				if(correctTest) data.err[e][i] = fabs(data.driver->innerproduct(data.mode[i],testmode));
				else            data.err[e][i] = nan("");
			}
			e++;
			if(!containsK) delete testmode;
		}
		//STEP 7c: for n=0 polytropes, we can compare frequencies to the exact Pekeris formula
		if(data.error[error::isIsopycnic]){
			for(int i=enext;i<data.mode_done; i++){
				if(data.k[i] >=0)
					data.err[e][i] = mode::compare_Pekeris(data.w[i], data.l[i], data.k[i], data.adiabatic_index);
				else data.err[e][i] = nan("");
			}
			e++;
		}
		//STEP 7d: for certain polytrope frequencies, we can compare to tables in JCD-DJM paper
		if(data.error[error::isJCD]) {
			double polytrope_index = data.input_params[0];
			for(int i=enext;i<data.mode_done; i++){
				if((data.l[i]==1 | data.l[i]==2 | data.l[i]==3) & (data.k[i]>0 & data.k[i]<36))
					data.err[e][i] = mode::compare_JCD(polytrope_index, data.l[i], data.k[i], data.w[i]);
				else data.err[e][i] = nan("");
			}
			e++;
		}
		if(e!=data.i_err) printf("Error in mode error listing...\n");
		printf("done\n");

//STEP 8: at the end of each L, print all data to the output file 
		io::write_mode_output(data);

		// clean up the remaining allocated pointers to modes
		// delete only the modes that were not in the list
		for(auto mt=modefilled.begin(); mt!=modefilled.end(); mt++){
			if(!kl.count(mt->second->modeOrder()))
				delete mt->second;
		}
	}

	// cleanup the sizes of L,K arrays to actually discovered data
	data.mode_num = nextmode;

	return 0;
}

} // namespace mode

#endif