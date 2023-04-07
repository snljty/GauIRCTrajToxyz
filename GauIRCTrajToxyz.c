# include <stdio.h>
# include <stdlib.h>
# include <string.h>

int const num_dims = 3;

void Print_help();

void Set_string_empty(char *s, int len);

void Set_buffer_string_empty(char *s);

void Check_command_arguments(int argc, char const **argv);

int Get_input_file_name(int argc, char const **argv, char *ifl_name_reverse, char *ifl_name_forward, char *ifl_name, char *line);

void Check_input_file_names(int num_input_files, char const *ifl_name_reverse, char const *ifl_name_forward, char const *ifl_name);

void Get_line_remove_end(char *line);

void Remove_line_end(char *line);

void Get_num_frames(int num_input_files, int *num_frames_reverse_ptr, int *num_frames_forward_ptr, \
    char const *ifl_name_reverse, char const *ifl_name_forward, char const *ifl_name, char *line);

void Get_num_atoms(int num_input_files, int *num_atoms_ptr, char const *ifl_name_reverse, char const *ifl_name_forward, char const *ifl_name, char *line);

double *Get_atom_coordinate_ptr(double *atom_coordinaes, int num_atoms, int ind_frame, int ind_atom);

inline double *Get_atom_coordinate_ptr(double *atom_coordinaes, int num_atoms, int ind_frame, int ind_atom)
{
    return atom_coordinaes + ind_frame * num_atoms * num_dims + ind_atom * num_dims;
}

/*
double Get_atom_coordinate_value(double const *atom_coordinaes, int num_atoms, int ind_frame, int ind_atom);

inline double Get_atom_coordinate_value(double const *atom_coordinaes, int num_atoms, int ind_frame, int ind_atom)
{
    return atom_coordinaes[ind_frame * num_atoms * num_dims + ind_atom * num_dims];
}
*/

void Get_atom_names_and_coordinates(int num_input_files, int num_frames_reverse, int num_frames_forward, int num_atoms, \
    char const *ifl_name_reverse, char const *ifl_name_forward, char const *ifl_name, char *line, \
    int *atom_indices, double *atom_coordinaes);

void Check_backup_name(char const *ofl_name, char const *ofl_name_bak);

void Write_traj(char const *ofl_name, int num_frames_reverse, int num_frames_forward, int num_atoms, \
    int const *atom_indices, double const *atom_coordinaes);

int main(int argc, char const **argv)
{
    /* ifl is used and only used if reverse and forward are in the same file */
    char ifl_name_reverse[BUFSIZ + 1] = "", ifl_name_forward[BUFSIZ + 1] = "", ifl_name[BUFSIZ + 1] = "";
    char line[BUFSIZ + 1] = "";
    int num_input_files;
    int num_atoms;
    int num_frames_reverse, num_frames_forward, num_frames;
    int *atom_indices = NULL;
    double *atom_coordinaes = NULL;
    char const ofl_name[] = "IRC_traj.xyz";
    char const ofl_name_bak[] = "IRC_traj.bak";

    num_input_files = Get_input_file_name(argc, argv, ifl_name_reverse, ifl_name_forward, ifl_name, line);
    Get_num_frames(num_input_files, & num_frames_reverse, & num_frames_forward, ifl_name_reverse, ifl_name_forward, ifl_name, line);
    num_frames = num_frames_reverse + 1 + num_frames_forward;
    Get_num_atoms(num_input_files, & num_atoms, ifl_name_reverse, ifl_name_forward, ifl_name, line);
    atom_indices = (int *)malloc(num_atoms * sizeof(int));
    atom_coordinaes = (double *)malloc(num_frames * num_atoms * num_dims * sizeof(double));
    Get_atom_names_and_coordinates(num_input_files, num_frames_reverse, num_frames_forward, num_atoms, \
        ifl_name_reverse, ifl_name_forward, ifl_name, line, atom_indices, atom_coordinaes);
    Check_backup_name(ofl_name, ofl_name_bak);
    Write_traj(ofl_name, num_frames_reverse, num_frames_forward, num_atoms, atom_indices, atom_coordinaes);

    free(atom_indices);
    atom_indices = NULL;
    free(atom_coordinaes);
    atom_coordinaes = NULL;

    return 0;
}

void Print_help(int argc, char const **argv)
{
    fprintf(stdout, "Usage: %s [reverse_IRC] [forward_IRC]\n", argv[0]);
    fprintf(stdout, "   or: %s [IRC]\n", argv[0]);
    fprintf(stdout, "Only supports Gaussian output file.\n");
    exit(EXIT_SUCCESS);

    return;
}

void Set_string_empty(char *s, int len)
{
    memset(s, 0, (len + 1));
    return;
}

void Set_buffer_string_empty(char *s)
{
    return Set_string_empty(s, BUFSIZ);
}

void Check_command_arguments(int argc, char const **argv)
{
    int iarg;

    if (argc - 1 > 2)
    {
        Print_help(argc, argv);
    }
    for (iarg = 1; iarg < argc; ++ iarg)
    {
        if (! (strcmp(argv[iarg], "--help") && strcmp(argv[iarg], "-h") && strcmp(argv[iarg], "/?")))
        {
            Print_help(argc, argv);
        }
    }

    return;
}

int Get_input_file_name(int argc, char const **argv, char *ifl_name_reverse, char *ifl_name_forward, char *ifl_name, char *line)
{
    /* returns the number of file names */
    int ret;

    Check_command_arguments(argc, argv);
    Set_buffer_string_empty(ifl_name_forward);
    Set_buffer_string_empty(ifl_name_reverse);
    Set_buffer_string_empty(ifl_name);
    if (! (argc - 1))
    {
        printf("Do you want to handle a IRC output file with the whole IRC path on both sized, \n");
        printf("or two seperated files with forward and reverse IRC path individually?\n");
        printf("1 for one file, 2 for two files.\n");
        Get_line_remove_end(line);
        if (sscanf(line, "%d", & ret) != 1 || ret != 1 && ret != 2)
        {
            fprintf(stderr, "Error! Cannot recognize \"%s\" as 1 or 2.\n", line);
            exit(EXIT_FAILURE);
        }
        if (ret == 1)
        {
            printf("Input the name of the IRC output file:\n");
            Get_line_remove_end(ifl_name);
        }
        else
        {
            printf("Input the name of the REVERSE IRC output file:\n");
            Get_line_remove_end(ifl_name_reverse);
            printf("Input the name of the FORWARD IRC output file:\n");
            Get_line_remove_end(ifl_name_forward);
        }
    }
    else
    {
        ret = argc - 1;
        if (ret == 1)
        {
            strncpy(ifl_name, argv[1], BUFSIZ + 1);
        }
        else
        {
            strncpy(ifl_name_reverse, argv[1], BUFSIZ + 1);
            strncpy(ifl_name_forward, argv[2], BUFSIZ + 1);
        }
    }
    Check_input_file_names(ret, ifl_name_reverse, ifl_name_forward, ifl_name);

    return ret;
}

void Check_input_file_names(int num_input_files, char const *ifl_name_reverse, char const *ifl_name_forward, char const *ifl_name)
{
    FILE *ifl = NULL;
    char *pos = NULL;

    if (num_input_files == 1)
    {
        pos = strrchr(ifl_name, '.');
        if (! pos)
        {
            fprintf(stderr, "Error! Cannot get the suffix of file name \"%s\".\n", ifl_name);
            exit(EXIT_FAILURE);
        }
        if (strcmp(pos, ".out") && strcmp(pos, ".log"))
        {
            fprintf(stderr, "Error! The suffix of the file name \"%s\" should be \".out\" or \".log\", but got \"%s\".\n", ifl_name, pos);
            exit(EXIT_FAILURE);
        }
        ifl = fopen(ifl_name, "rt");
        if (! ifl)
        {
            fprintf(stderr, "Error! Cannot open \"%s\" for reading.\n", ifl_name);
            exit(EXIT_FAILURE);
        }
        fclose(ifl);
        ifl = NULL;
    }
    else
    {
        pos = strrchr(ifl_name_reverse, '.');
        if (! pos)
        {
            fprintf(stderr, "Error! Cannot get the suffix of file name \"%s\".\n", ifl_name_reverse);
            exit(EXIT_FAILURE);
        }
        if (strcmp(pos, ".out") && strcmp(pos, ".log"))
        {
            fprintf(stderr, "Error! The suffix of the file name \"%s\" should be \".out\" or \".log\", but got \"%s\".\n", ifl_name_reverse, pos);
            exit(EXIT_FAILURE);
        }
        pos = strrchr(ifl_name_forward, '.');
        if (! pos)
        {
            fprintf(stderr, "Error! Cannot get the suffix of file name \"%s\".\n", ifl_name_forward);
            exit(EXIT_FAILURE);
        }
        if (strcmp(pos, ".out") && strcmp(pos, ".log"))
        {
            fprintf(stderr, "Error! The suffix of the file name \"%s\" should be \".out\" or \".log\", but got \"%s\".\n", ifl_name_forward, pos);
            exit(EXIT_FAILURE);
        }
        ifl = fopen(ifl_name_reverse, "rt");
        if (! ifl)
        {
            fprintf(stderr, "Error! Cannot open \"%s\" for reading.\n", ifl_name_reverse);
            exit(EXIT_FAILURE);
        }
        fclose(ifl);
        ifl = NULL;
        ifl = fopen(ifl_name_forward, "rt");
        if (! ifl)
        {
            fprintf(stderr, "Error! Cannot open \"%s\" for reading.\n", ifl_name_forward);
            exit(EXIT_FAILURE);
        }
        fclose(ifl);
        ifl = NULL;
    }

    return;
}

void Get_line_remove_end(char *line)
{
    if (! fgets(line, BUFSIZ, stdin))
    {
        fprintf(stderr, "Error! Cannot get a line as input.\n");
        exit(EXIT_FAILURE);
    }
    return Remove_line_end(line);
}

void Remove_line_end(char *line)
{
    int len = strlen(line);

    if (len)
    {
        if (line[len - 1] == '\n')
        {
            line[len - 1] = '\0';
        }
    }

    return;
}

void Get_num_frames(int num_input_files, int *num_frames_reverse_ptr, int *num_frames_forward_ptr, \
    char const *ifl_name_reverse, char const *ifl_name_forward, char const *ifl_name, char *line)
{
    FILE *ifl = NULL;

    if (num_input_files == 1)
    {
        * num_frames_forward_ptr = -1;
        ifl = fopen(ifl_name, "rt");
        while (fgets(line, BUFSIZ, ifl))
        {
            if (strstr(line, "Calculation of FORWARD path complete."))
            {
                break;
            }
            if (strstr(line, "Path Number:   1"))
            {
                ++ * num_frames_forward_ptr;
            }
        }
        * num_frames_reverse_ptr = 0;
        while (fgets(line, BUFSIZ, ifl))
        {
            if (strstr(line, "Path Number:   2"))
            {
                ++ * num_frames_reverse_ptr;
            }
        }
        fclose(ifl);
        ifl = NULL;
    }
    else
    {
        * num_frames_reverse_ptr = -1;
        ifl = fopen(ifl_name_reverse, "rt");
        while (fgets(line, BUFSIZ, ifl))
        {
            if (strstr(line, "Path Number:   1"))
            {
                ++ * num_frames_reverse_ptr;
            }
        }
        fclose(ifl);
        ifl = NULL;
        * num_frames_forward_ptr = -1;
        ifl = fopen(ifl_name_forward, "rt");
        while (fgets(line, BUFSIZ, ifl))
        {
            if (strstr(line, "Path Number:   1"))
            {
                ++ * num_frames_forward_ptr;
            }
        }
        fclose(ifl);
        ifl = NULL;
    }

    return;
}

void Get_num_atoms(int num_input_files, int *num_atoms_ptr, char const *ifl_name_reverse, char const *ifl_name_forward, char const *ifl_name, char *line)
{
    FILE *ifl = NULL;
    char *pos = NULL;
    int i;

    ifl = num_input_files == 1 ? fopen(ifl_name, "rt") : fopen(ifl_name_reverse, "rt");
    while (fgets(line, BUFSIZ, ifl))
    {
        pos = strstr(line, "NAtoms=");
        if (pos)
        {
            sscanf(pos + strlen("NAtoms="), "%d", num_atoms_ptr);
            break;
        }
    }
    fclose(ifl);
    ifl = NULL;

    return;
}

void Get_atom_names_and_coordinates(int num_input_files, int num_frames_reverse, int num_frames_forward, int num_atoms, \
    char const *ifl_name_reverse, char const *ifl_name_forward, char const *ifl_name, char *line, \
    int *atom_indices, double *atom_coordinaes)
{
    int ind_frame, ind_atom;
    FILE *ifl = NULL;
    int i;
    int const num_useless_lines = 4;
    double *coord_pos = NULL;

    /* read the last "Input orientation:" before "Point Number:", WTF */
    if (num_input_files == 1)
    {
        ifl = fopen(ifl_name, "rt");
        /* TS and forward */
        for (ind_frame = num_frames_reverse; ind_frame < num_frames_reverse + 1 + num_frames_forward; ++ ind_frame)
        {
            while (fgets(line, BUFSIZ, ifl))
            {
                if (strstr(line, "Path Number:   1"))
                {
                    break;
                }
                if (strstr(line, "Input orientation:"))
                {
                    for (i = 0; i < num_useless_lines; ++ i)
                    {
                        fgets(line, BUFSIZ, ifl);
                    }
                    for (ind_atom = 0; ind_atom < num_atoms; ++ ind_atom)
                    {
                        fgets(line, BUFSIZ, ifl);
                        coord_pos = Get_atom_coordinate_ptr(atom_coordinaes, num_atoms, ind_frame, ind_atom);
                        sscanf(line, "%*d%d%*d%lg%lg%lg", atom_indices + ind_atom, coord_pos, coord_pos + 1, coord_pos + 2);
                    }
                }
            }
        }
        /* reverse */
        for (ind_frame = num_frames_reverse - 1; ind_frame >= 0; -- ind_frame)
        {
            while (fgets(line, BUFSIZ, ifl))
            {
                if (strstr(line, "Path Number:   2"))
                {
                    break;
                }
                if (strstr(line, "Input orientation:"))
                {
                    for (i = 0; i < num_useless_lines; ++ i)
                    {
                        fgets(line, BUFSIZ, ifl);
                    }
                    for (ind_atom = 0; ind_atom < num_atoms; ++ ind_atom)
                    {
                        fgets(line, BUFSIZ, ifl);
                        coord_pos = Get_atom_coordinate_ptr(atom_coordinaes, num_atoms, ind_frame, ind_atom);
                        sscanf(line, "%*d%d%*d%lg%lg%lg", atom_indices + ind_atom, coord_pos, coord_pos + 1, coord_pos + 2);
                    }
                }
            }
        }

        fclose(ifl);
        ifl = NULL;
    }
    else
    {
        ifl = fopen(ifl_name_reverse, "rt");
        for (ind_frame = num_frames_reverse; ind_frame >= 0; -- ind_frame)
        {
            while (fgets(line, BUFSIZ, ifl))
            {
                if (strstr(line, "Path Number:   1"))
                {
                    break;
                }
                if (strstr(line, "Input orientation:"))
                {
                    for (i = 0; i < num_useless_lines; ++ i)
                    {
                        fgets(line, BUFSIZ, ifl);
                    }
                    for (ind_atom = 0; ind_atom < num_atoms; ++ ind_atom)
                    {
                        fgets(line, BUFSIZ, ifl);
                        coord_pos = Get_atom_coordinate_ptr(atom_coordinaes, num_atoms, ind_frame, ind_atom);
                        sscanf(line, "%*d%d%*d%lg%lg%lg", atom_indices + ind_atom, coord_pos, coord_pos + 1, coord_pos + 2);
                    }
                }
            }
        }
        fclose(ifl);
        ifl = NULL;
        ifl = fopen(ifl_name_forward, "rt");
        for (ind_frame = num_frames_reverse; ind_frame < num_frames_reverse + 1 + num_frames_forward; ++ ind_frame)
        {
            while (fgets(line, BUFSIZ, ifl))
            {
                if (strstr(line, "Path Number:   1"))
                {
                    break;
                }
                if (strstr(line, "Input orientation:"))
                {
                    for (i = 0; i < num_useless_lines; ++ i)
                    {
                        fgets(line, BUFSIZ, ifl);
                    }
                    for (ind_atom = 0; ind_atom < num_atoms; ++ ind_atom)
                    {
                        fgets(line, BUFSIZ, ifl);
                        coord_pos = Get_atom_coordinate_ptr(atom_coordinaes, num_atoms, ind_frame, ind_atom);
                        sscanf(line, "%*d%d%*d%lg%lg%lg", atom_indices + ind_atom, coord_pos, coord_pos + 1, coord_pos + 2);
                    }
                }
            }
        }
        fclose(ifl);
        ifl = NULL;
    }

    return;
}

void Check_backup_name(char const *ofl_name, char const *ofl_name_bak)
{
    FILE *ofl = NULL;

    ofl = fopen(ofl_name, "rt");
    if (ofl)
    {
        fclose(ofl);
        ofl = NULL;
        fprintf(stderr, "Warning: file \"%s\" already exists, will be moved to \"%s\".\n", ofl_name, ofl_name_bak);
        ofl = fopen(ofl_name_bak, "rt");
        if (ofl)
        {
            fclose(ofl);
            ofl = NULL;
            fprintf(stderr, "Warning: backup file \"%s\" already exists, will be overwritten.\n", ofl_name_bak);
            remove(ofl_name_bak);
        }
        rename(ofl_name, ofl_name_bak);
    }

    return;
}

void Write_traj(char const *ofl_name, int num_frames_reverse, int num_frames_forward, int num_atoms, \
    int const *atom_indices, double const *atom_coordinaes)
{
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
    int ind_frame, ind_atom;
    FILE *ofl = NULL;
    double const* coord_pos = NULL;

    ofl = fopen(ofl_name, "wt");
    for (ind_frame = 0; ind_frame < num_frames_reverse; ++ ind_frame)
    {
        fprintf(ofl, "%d\n", num_atoms);
        fprintf(ofl, "TS %+d\n", ind_frame - num_frames_reverse);
        for (ind_atom = 0; ind_atom < num_atoms; ++ ind_atom)
        {
            coord_pos = Get_atom_coordinate_ptr((double *)atom_coordinaes, num_atoms, ind_frame, ind_atom);
            fprintf(ofl, " %-2s        %10.6lf    %10.6lf    %10.6lf\n", elements_list[atom_indices[ind_atom]], \
                coord_pos[0], coord_pos[1], coord_pos[2]);
        }
    }
    ind_frame = num_frames_reverse;
    fprintf(ofl, "%d\n", num_atoms);
    fprintf(ofl, "TS\n");
    for (ind_atom = 0; ind_atom < num_atoms; ++ ind_atom)
    {
        coord_pos = Get_atom_coordinate_ptr((double *)atom_coordinaes, num_atoms, ind_frame, ind_atom);
        fprintf(ofl, " %-2s        %10.6lf    %10.6lf    %10.6lf\n", elements_list[atom_indices[ind_atom]], \
            coord_pos[0], coord_pos[1], coord_pos[2]);
    }
    for (ind_frame = num_frames_reverse + 1; ind_frame < num_frames_reverse + 1 + num_frames_forward; ++ ind_frame)
    {
        fprintf(ofl, "%d\n", num_atoms);
        fprintf(ofl, "TS %+d\n", ind_frame - num_frames_reverse);
        for (ind_atom = 0; ind_atom < num_atoms; ++ ind_atom)
        {
            coord_pos = Get_atom_coordinate_ptr((double *)atom_coordinaes, num_atoms, ind_frame, ind_atom);
            fprintf(ofl, " %-2s        %10.6lf    %10.6lf    %10.6lf\n", elements_list[atom_indices[ind_atom]], \
                coord_pos[0], coord_pos[1], coord_pos[2]);
        }
    }
    fclose(ofl);
    ofl = NULL;
    printf("Done. Written to \"%s\".\n", ofl_name);

    return;
}

