/* works only for Gaussian 16 */

/* frame 1 is TS, then until "Calculation of FORWARD path complete." is forward, */
/* from "Beginning calculation of the REVERSE path." to "Calculation of REVERSE path complete." is reverse. */

/* There is an extra coordinates output after the complection of the IRC, which is the same as the 1st frame of reverse */

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdbool.h>

int main(int argc, char const *argv[])
{
    # define max_ene_type_str_len 10
    # define max_coordinates_locator_len 25
    # define num_useless_lines_after_coordinates_locator 4
    char const *elements_list[] = {"", \
     "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , \
     "F" , "Ne", "Na", "Mg", "Al", "Si", "P" , "S" , \
     "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", \
     "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", \
     "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", \
     "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", \
     "In", "Sn", "Sb", "Te", "I" , "Xe", "Cs", "Ba", \
     "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", \
     "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", \
     "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg", \
     "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", \
     "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am", "Cm", \
     "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", \
     "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", \
     "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
    char buf[BUFSIZ + 1] = "";
    size_t *read_pos = NULL;
    int pos_len_allocated = 0;
    char *tok = NULL;
    char ifl_name[BUFSIZ + 1] = "";
    char const ofl_name[] = "IRCsplit.xyz";
    char const ofl_bak_name[] = "IRCsplit.bak";
    FILE *ifl = NULL;
    FILE *ofl = NULL;
    char coordinates_locator[max_coordinates_locator_len + 1] = "";
    char ene_type_str[max_ene_type_str_len + 1] = "";
    double ene = 0.0;
    int num_frames = 0;
    int ind_frame_read = 0;
    int ind_frame_write = 0;
    int ind_reverse_beg_frame = 0;
    double x = 0.0, y = 0.0, z = 0.0;
    char c = '\0';
    int element_index = 0;
    int num_atoms = 0;
    int i = 0;
    bool is_ts_reached = false;

    /* initialize variables */
    tok = NULL;
    ifl = NULL;
    ofl = NULL;
    read_pos = NULL;
    memset(buf, 0, (BUFSIZ + 1) * sizeof(char));
    memset(ifl_name, 0, (BUFSIZ + 1) * sizeof(char));
    memset(coordinates_locator, 0, (max_coordinates_locator_len + 1) * sizeof(char));
    memset(ene_type_str, 0, (max_ene_type_str_len + 1) * sizeof(char));
    ene = 0.0;
    num_frames = 0;
    ind_frame_read = 0;
    ind_frame_write = 0;
    ind_reverse_beg_frame = 0;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    c = '\0';
    element_index = 0;
    num_atoms = 0;
    i = 0;
    is_ts_reached = false;

    /* pharse command arguments, get file names */
    if (argc - 1 >= 2)
    {
        fprintf(stderr, "Error! Too many command arguments! At most 1 required, but got %d.\n", argc - 1);
        exit(EXIT_FAILURE);
    }
    else if (argc - 1 == 1)
    {
        strncpy(ifl_name, argv[1], BUFSIZ + 1);
    }
    else /* argc - 1 == 0 */
    {
        printf("Input file name of output of a Gaussian IRC task: \n");
        while (! fgets(buf, BUFSIZ, stdin))
        {
            ;
        }
        if (buf[strlen(buf) - 1] == '\n')
        {
            buf[strlen(buf) - 1] = '\0';
        }
        if (buf[0] == '\"')
        {
            buf[strlen(buf) - 1] = '\0';
            strcpy(ifl_name, buf + 1);
        }
        else
        {
            strcpy(ifl_name, buf);
        }
    }

    /* check suffix of input file name */
    if (! strrchr(ifl_name, '.') || strcmp(strrchr(ifl_name, '.'), ".out") && strcmp(strrchr(ifl_name, '.'), ".log"))
    {
        fprintf(stderr, "Error! The suffix of the input file name should be either \".out\" or \".log\", but it is not.\n");
        exit(EXIT_FAILURE);
    }

    /* check input and output files, and open them */
    if (! (ifl = fopen(ifl_name, "rt")))
    {
        fprintf(stderr, "Error! Cannot read \"%s\".\n", ifl_name);
        exit(EXIT_FAILURE);
    }
    if (ofl = fopen(ofl_name, "rt"))
    {
        fprintf(stderr, "Warning! File \"%s\" already exists. Will try to move it to \"%s\".\n", ofl_name, ofl_bak_name);
        fclose(ofl);
        ofl = NULL;
        if (ofl = fopen(ofl_bak_name, "rt"))
        {
            fprintf(stderr, "Warning! Old backup file \"%s\" already exists, and it will be removed.\n", ofl_bak_name);
            fclose(ofl);
            ofl = NULL;
            remove(ofl_bak_name);
        }
        rename(ofl_name, ofl_bak_name);
    }
    ofl = fopen(ofl_name, "wt");

    /* get amount of atoms */
    while (fgets(buf, BUFSIZ, ifl))
    {
        if (tok = strstr(buf, "NAtoms="))
        {
            tok += strlen("NAtoms=");
            sscanf(tok, "%d", & num_atoms);
        }
    }
    rewind(ifl);

    /* get coordinate locator */
    while (true)
    {
        if (! fgets(buf, BUFSIZ, ifl))
        {
            fprintf(stderr, "Error! Cannot find \"Input\", \"Standard\" or \"Z-Matrix\" orientation.\n");
            fclose(ifl);
            ifl = NULL;
            fclose(ofl);
            ofl = NULL;
            exit(EXIT_FAILURE);
        }
        if (strstr(buf, "Input orientation"))
        {
            strcpy(coordinates_locator, "Input orientation");
            break;
        }
        if (strstr(buf, "Standard orientation"))
        {
            strcpy(coordinates_locator, "Standard orientation");
            break;
        }
        if (strstr(buf, "Z-Matrix orientation"))
        {
            strcpy(coordinates_locator, "Z-Matrix orientation");
            break;
        }
    }
    rewind(ifl);

    /* check energy type */
    while (true)
    {
        if (! fgets(buf, BUFSIZ, ifl))
        {
            fprintf(stderr, "Error! Cannot determine energy type.\n");
            fclose(ifl);
            ifl = NULL;
            fclose(ofl);
            ofl = NULL;
            exit(EXIT_FAILURE);
        }
        if (strstr(buf, "Energy="))
        {
            strcpy(ene_type_str, "MM");
            break;
        }
        if (strstr(buf, "SCF Done"))
        {
            break;
        }
    }
    if (strcmp(ene_type_str, "MM")) /* found "SCF Done" */
    {
        while (true)
        {
            if (! fgets(buf, BUFSIZ, ifl) || strstr(buf, "Population analysis"))
            {
                strcpy(ene_type_str, "SCF");
                break;
            }
            if (strstr(buf, "EUMP2 ="))
            {
                strcpy(ene_type_str, "MP2");
                break;
            }
            if (! strncmp(buf, " E2(", strlen(" E2(")))
            {
                strcpy(ene_type_str, "DFTPT2");
                break;
            }
            if (strstr(buf, "E(CIS/TDA)"))
            {
                strcpy(ene_type_str, "CIS/TDA");
                break;
            }
            if (strstr(buf, "E(TD-HF/TD-DFT)"))
            {
                strcpy(ene_type_str, "TD");
                break;
            }
        }
    }
    rewind(ifl);

    /* get amount of frames */
    num_frames = 0;
    while (fgets(buf, BUFSIZ, ifl))
    {
        /* if (strstr(buf, "Summary of reaction path following")) */
        if (strstr(buf, "Calculation of REVERSE path complete."))
        {
            break;
        }
        if (strstr(buf, "Calculation of FORWARD path complete."))
        {
            ind_reverse_beg_frame = num_frames + 1;
        }
        if (strstr(buf, coordinates_locator))
        {
            ++ num_frames;
        }
    }
    rewind(ifl);

    /* get positions of each frames */
    /* read_pos[ind_frame_read - 1] is for frame ind_frame_read */
    read_pos = (size_t *)malloc(num_frames * sizeof(size_t));
    ind_frame_read = 0;
    while (fgets(buf, BUFSIZ, ifl))
    {
        /* if (strstr(buf, "Summary of reaction path following")) */
        if (strstr(buf, "Calculation of REVERSE path complete."))
        {
            break;
        }
        if (strstr(buf, coordinates_locator))
        {
            ++ ind_frame_read;
            read_pos[ind_frame_read - 1] = ftell(ifl);
        }
    }
    rewind(ifl);

    /* read and write */
    ind_frame_read = num_frames;
    is_ts_reached = false;
    for (ind_frame_write = 1; ind_frame_write <= num_frames; ++ ind_frame_write)
    {
        fseek(ifl, read_pos[ind_frame_read - 1], SEEK_SET);
        /* read energy here */
        if (strstr(ene_type_str, "MM"))
        {
            while (fgets(buf, BUFSIZ, ifl))
            {
                if (tok = strstr(buf, "Energy="))
                {
                    break;
                }
            }
            tok += strlen("Energy=");
            sscanf(tok, "%lg", & ene);
        }
        else if (strstr(ene_type_str, "SCF"))
        {
            while (fgets(buf, BUFSIZ, ifl))
            {
                if (strstr(buf, "SCF Done"))
                {
                    break;
                }
            }
            tok = strchr(buf, '=') + strlen("=");
            sscanf(tok, "%lg", & ene);
        }
        else if (strstr(ene_type_str, "MP2"))
        {
            while (fgets(buf, BUFSIZ, ifl))
            {
                if (tok = strstr(buf, "EUMP2"))
                {
                    break;
                }
            }
            tok = strchr(tok, ' ') + strlen("=");
            * strchr(tok, 'D') = 'E';
            sscanf(tok, "%lg", & ene);
        }
        else if (strstr(ene_type_str, "DFTPT2"))
        {
            while (fgets(buf, BUFSIZ, ifl))
            {
                if (! strncmp(buf, " E2(", strlen(" E2(")))
                {
                    break;
                }
            }
            tok = strstr(buf, "E(");
            tok = strchr(tok, '=') + strlen("=");
            * strchr(tok, 'D') = 'E';
            sscanf(tok, "%lg", & ene);
        }
        else if (strstr(ene_type_str, "CIS/TDA"))
        {
            while (fgets(buf, BUFSIZ, ifl))
            {
                if (tok = strstr(buf, "E(CIS/TDA)"))
                {
                    break;
                }
            }
            tok = strchr(buf, '=') + strlen("=");
            sscanf(tok, "%lg", & ene);
        }
        else if (strstr(ene_type_str, "TD"))
        {
            while (fgets(buf, BUFSIZ, ifl))
            {
                if (tok = strstr(buf, "E(TD-HF/TD-DFT)"))
                {
                    break;
                }
            }
            tok = strchr(buf, '=') + strlen("=");
            sscanf(tok, "%lg", & ene);
        }
        else
        {
            /* should never happen */
            ;
        }
        tok = NULL;
        fseek(ifl, read_pos[ind_frame_read - 1], SEEK_SET);
        /* read geometry here */
        fprintf(ofl, "%4d\n", num_atoms);
        if (ind_frame_read >= ind_reverse_beg_frame)
        {
            fprintf(ofl, "frame %4d (<-): energy = %17.10lf Hartree (\"%s\" level)\n", ind_frame_write, ene, ene_type_str);
        }
        else if (ind_frame_read > 1)
        {
            fprintf(ofl, "frame %4d (->): energy = %17.10lf Hartree (\"%s\" level)\n", ind_frame_write, ene, ene_type_str);
        }
        else
        {
            fprintf(ofl, "frame %4d (TS): energy = %17.10lf Hartree (\"%s\" level)\n", ind_frame_write, ene, ene_type_str);
        }
        for (i = 0; i < num_useless_lines_after_coordinates_locator; ++ i)
        {
            fgets(buf, BUFSIZ, ifl);
        }
        for (i = 0; i < num_atoms; ++ i)
        {
            fgets(buf, BUFSIZ, ifl);
            sscanf(buf, "%*d%d%*d%lg%lg%lg", & element_index, & x, & y, & z);
            fprintf(ofl, " %-2s %15s %12.8lf    %12.8lf    %12.8lf\n", elements_list[element_index], "", x, y, z);
        }
        fflush(ofl);
        if (ind_frame_read == ind_reverse_beg_frame)
        {
            ind_frame_read = 1;
            /* reverse from num_frames to ind_reverse_beg_frame, then 1, then forward from */
            /* 2 to ind_reverse_beg_frame - 1. ind_reverse_beg_frame - 1 == ind_forward_end_frame */
        }
        else
        {
            ind_frame_read += ind_frame_read > ind_reverse_beg_frame ? -1 : 1;
        }
    }

    /* close all files and release memory */
    free(read_pos);
    read_pos = NULL;
    fclose(ifl);
    ifl = NULL;
    fclose(ofl);
    ofl = NULL;

    /* final output */
    printf("Input file name: \"%s\".\n", ifl_name);
    printf("Number of atoms: %d.\n", num_atoms);
    printf("Number of frames: %d.\n", num_frames);
    printf("Index of the frame correspond to transition state: %d.\n", num_frames - ind_reverse_beg_frame + 2);
    printf("The result has been saved to \"%s\".\n", ofl_name);
    printf("Everything is done!\n");
    if (! (argc - 1))
    {
        printf("Press <Enter> to exit ...");
        while ((c = getchar()) != '\n' && c != EOF)
        {
            ;
        }
    }

    return 0;
}

