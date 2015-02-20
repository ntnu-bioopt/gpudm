#ifndef CHROMOPHORES_H_DEFINED
#define CHROMOPHORES_H_DEFINED
#include <string>

class Chromophores{
	public:
		Chromophores(); //blood is set as standard
		void checkAndSet(bool *val);
		void checkAndUnset(bool *val);

		void setWat(){checkAndSet(&haswat);};
		void setKonst(){checkAndSet(&haskonst);};
		void setMel(){checkAndSet(&hasmel);};
		void setMethb(){checkAndSet(&hasmethb);};
		void setBil(){checkAndSet(&hasbil);};
		void setBet(){checkAndSet(&hasbet);};
		void setDeoxy(){checkAndSet(&hasdeoxy);};
		void setOxy(){checkAndSet(&hasoxy);};
		void setCOHb(){checkAndSet(&hascohb);};
		void setFat(){checkAndSet(&hasfat);};

		void unsetWat(){checkAndUnset(&haswat);};
		void unsetKonst(){checkAndUnset(&haskonst);};
		void unsetMel(){checkAndUnset(&hasmel);};
		void unsetMethb(){checkAndUnset(&hasmethb);};
		void unsetBil(){checkAndUnset(&hasbil);};
		void unsetBet(){checkAndUnset(&hasbet);};
		void unsetDeoxy(){checkAndUnset(&hasdeoxy);};
		void unsetOxy(){checkAndUnset(&hasoxy);};
		void unsetCOHb(){checkAndUnset(&hascohb);};
		void unsetFat(){checkAndUnset(&hasfat);};

		int getNumEndmembers(){return numEnd;};
		float *getAbsArray(float w);
		float getAbs(float w, int end);
		std::string getName(int end);
		int getMelInd();
		int getMethbInd();
		int getKonstInd();


		float getMaxVal(int endmember);
		float getMinVal(int endmember);

		bool hasOxy(){return hasoxy;};
		bool hasDeoxy(){return hasdeoxy;};
	private:
		bool hasoxy;
		bool hasdeoxy;
		bool haswat;
		bool haskonst;
		bool hasmel;
		bool hasmethb;
		bool hasbil;
		bool hasbet;
		bool hascohb;
		bool hasfat;

		int numEnd;
		
		float oxymax, deoxymax, methbmax, bilmax, betmax, watmax, cohbmax, melmax, konstmax, fatmax;
};

#endif
