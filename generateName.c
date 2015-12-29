/*Auxillary function file to generate the corresponding filenames.
	case 1-4 generates 1 filename which includes the following data:
		
		(i)The External Magnetic field Hex vector
		(ii)Initial Values of theta and phi
		(iii)Value of current if necessary in case of Spin Transfer Torque and Field Like torque (cases 3 & 4)


	case 5-8 creates 2 filenames
		The first filename contains the the following data:
		(i) and (ii) are same as for above
		(iii) Value of current in case of Spin Transfer Torque and Field Like Torque (cases 6 & 7)
		(iv) The number of Ensembles used

		The second filename is for the file used to store the switching times of each ensemble.
		It has the following data:
		(i) and (ii) are the same as above.
		(iii) Value of current used.
		(iv) No of Ensembles
*/



char **generateFileNames(int op)
{

	
	char **filename;
	switch(op)
	{
		case 1:
		{
			filename = malloc(1*sizeof(char*));
			*filename = malloc(150*sizeof(char));
			strcat(*filename, "hexonly_Hex=");
			char *c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*filename, c);
			free(c);
			strcat(*filename, "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			/*c= malloc(20*sizeof(char));
			strcat(*filename, "_current=")
			sprintf(c, "%e", current);
			strcat(*filename, c);
			free(c);*/
			strcat(*(filename+1), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename+1), c);
			free(c);
			strcat(*filename, ".d");
			break;
		}

		case 2:
		{
			filename = malloc(1*sizeof(char*));
			*filename = malloc(150*sizeof(char));
			strcat(*filename, "hexDamp_Hex=");
			char *c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*filename, c);
			free(c);
			strcat(*filename, "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			strcat(*(filename+1), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename+1), c);
			free(c);
			strcat(*filename, ".d");
			break;
		}

		case 3:
		{
			filename = malloc(1*sizeof(char*));
			*filename = malloc(200*sizeof(char));
			strcat(*filename, "hexDampStt_Hex=");
			char *c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*filename, c);
			free(c);
			strcat(*filename, "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*filename, "_current=");
			sprintf(c, "%e", current);
			strcat(*filename, c);
			free(c);
			strcat(*(filename+1), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename+1), c);
			free(c);
			strcat(*filename, ".d");
			break;
		}

		case 4:
		{
			filename = malloc(1*sizeof(char*));
			*filename = malloc(250*sizeof(char));
			strcat(*filename, "hexDampSttFlt_Hex=");
			char *c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*filename, c);
			free(c);

			strcat(*filename, "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*filename, "_current=");
			sprintf(c, "%e", current);
			strcat(*filename, c);
			free(c);
			strcat(*(filename), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename), c);
			free(c);
			printf("Good till here\n");
			strcat(*filename, ".d");
			break;
		}

		case 5:
		{

			filename = malloc(2*sizeof(char*));
			*filename = malloc(250*sizeof(char));
			strcat(*filename, "hexDampFluc_Hex=");
			char *c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*filename, c);
			free(c);
			strcat(*filename, "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*filename, "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*filename, c);
			free(c);
			/*c= malloc(20*sizeof(char));
			strcat(*filename, "_current=")
			sprintf(c, "%e", current);
			strcat(*filename, c);
			free(c);*/
			strcat(*(filename), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename), c);
			free(c);
			strcat(*filename, ".d");


			*(filename+1) = malloc(250*sizeof(char));
			strcat(*(filename+1), "probDatahexDampFluc");
			c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*(filename+1), c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*(filename+1), "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*(filename+1), c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*(filename+1), "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), ".d");
	
			/*&*(filename+2) = malloc(150*sizeof(char));
			strcat(*(filename+2), "dwellTimeDatahexDampFluc");
			c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*(filename+2), c);
			free(c);
			strcat(*(filename+2), "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*(filename+2), c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*(filename+2), "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*(filename+2), c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*(filename+2), "_Ensembles=");
			sprintf(c, "%d", NO_OF_ENSEMBLE);
			strcat(*(filename+2), c);
			free(c);
			strcat(*(filename+2), ".d");	*/		
			break;
		}

		case 6:
		{
			filename = malloc(2*sizeof(char*));
			*filename = malloc(250*sizeof(char));
			strcat(*filename, "hexDampSttFluc_Hex=");
			char *c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*filename, c);
			free(c);
			strcat(*filename, "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*filename, "_current=");
			sprintf(c, "%e", current);
			strcat(*filename, c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*filename, "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*filename, c);
			free(c);
			strcat(*(filename), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename), c);
			free(c);
			strcat(*filename, ".d");


			*(filename+1) = malloc(250*sizeof(char));
			strcat(*(filename+1), "probDatahexDampSttFluc");
			c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*(filename+1), c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*(filename+1), "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*(filename+1), c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*(filename+1), "_current=");
			sprintf(c, "%e", current);
			strcat(*(filename+1), c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*(filename+1), "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), ".d");

			
			/**(filename+2) = malloc(150*sizeof(char));
			strcat(*(filename+2), "dwellTimeDatahexDampSttFluc");
			c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*(filename+2), c);
			free(c);
			strcat(*(filename+2), "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*(filename+2), c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*(filename+2), "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*(filename+2), c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*(filename+2), "_Ensembles=");
			sprintf(c, "%d", NO_OF_ENSEMBLE);
			strcat(*(filename+2), c);
			free(c);
			strcat(*(filename+2), ".d");*/			
			break;
		}

		case 7:
		{
			
			filename = malloc(2*sizeof(char*));
			*filename = malloc(250*sizeof(char));
			strcat(*filename, "hexDampSttFltFluc_Hex=");
			char *c = malloc(20*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*filename, c);
			free(c);
			
			strcat(*filename, "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*filename, "_current=");
			sprintf(c, "%e", current);
			strcat(*filename, c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*filename, "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*filename, c);
			free(c);
			strcat(*filename, "_Mf=");			
			c = malloc(25*sizeof(char));
			//printf("creating files\n");
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*filename, c);
			//printf("creating files2\n");
			free(c);
			strcat(*filename, ".d");



			*(filename+1) = malloc(250*sizeof(char));
			strcat(*(filename+1), "probDatahexDampSttFltFluc");
			c = malloc(20*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*(filename+1), c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*(filename+1), "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*(filename+1), c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*(filename+1), "_current=");
			sprintf(c, "%e", current);
			strcat(*(filename+1), c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*(filename+1), "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), "_Mf=");			
			c = malloc(25*sizeof(char));
			sprintf(c, "(%6.2f, %6.2f, %6.2f)", (Mf[0]), (Mf[1]), (Mf[2]));
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), ".d");
	
			
			/**(filename+2) = malloc(150*sizeof(char));
			strcat(*(filename+2), "dwellTimeDatahexDampSttFltFluc");
			c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*(filename+2), c);
			free(c);
			strcat(*(filename+2), "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*(filename+2), c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*(filename+2), "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*(filename+2), c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*(filename+2), "_Ensembles=");
			sprintf(c, "%d", NO_OF_ENSEMBLE);
			strcat(*(filename+2), c);
			free(c);
			strcat(*(filename+2), ".d");*/			
			break;
		}

		case 8:
		{
			filename = malloc(1*sizeof(char*));
			*filename = malloc(150*sizeof(char));
			strcat(*filename, "flucOnly");
			strcat(*filename, "_theta_initial=");
			char *c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*filename, "_current=");
			sprintf(c, "%e", current);
			strcat(*filename, c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*filename, "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*filename, c);
			free(c);
			strcat(*filename, ".d");			
			break;
		}

		default:
		{
			filename = malloc(3*sizeof(char*));
			*filename = malloc(150*sizeof(char));
			strcat(*filename, "hexDampSttFltFluc_Hex=");
			char *c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*filename, c);
			free(c);
			strcat(*filename, "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*filename, c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*filename, "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*filename, c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*filename, "_current=");
			sprintf(c, "%e", current);
			strcat(*filename, c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*filename, "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*filename, c);
			free(c);
			strcat(*filename, ".d");


			*(filename+1) = malloc(100*sizeof(char));
			strcat(*(filename+1), "probDatahexDampSttFltFluc");
			c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*(filename+1), c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*(filename+1), "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*(filename+1), c);
			free(c);
			c= malloc(20*sizeof(char));
			strcat(*filename, "_current=");
			sprintf(c, "%e", current);
			strcat(*(filename+1), c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*(filename+1), "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*(filename+1), c);
			free(c);
			strcat(*(filename+1), ".d");	
	
			
			*(filename+2) = malloc(150*sizeof(char));
			strcat(*(filename+2), "dwellTimeDatahexDampSttFltFluc");
			c = malloc(15*sizeof(char));
			sprintf(c, "(%d, %d, %d)", (int)(hex[0]*4*M_PI*Ms), (int)(hex[1]*4*M_PI*Ms), (int)(hex[2]*4*M_PI*Ms));
			strcat(*(filename+2), c);
			free(c);
			strcat(*(filename+2), "_theta_initial=");
			c = malloc(10*sizeof(char));
			sprintf(c, "%.4f", theta);
			strcat(*(filename+2), c);
			free(c);
			c = malloc(10*sizeof(char));
			strcat(*(filename+2), "_phi_initial=");
			sprintf(c, "%.4f", phi);
			strcat(*(filename+2), c);
			free(c);
			c = malloc(15*sizeof(char));
			strcat(*(filename+2), "_Ensembles=");
			sprintf(c, "%ld", NO_OF_ENSEMBLE);
			strcat(*(filename+2), c);
			free(c);
			strcat(*(filename+2), ".d");		
			break;
		}	
		
	}

	return filename;

}
