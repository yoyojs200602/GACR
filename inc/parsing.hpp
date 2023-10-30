#pragma once

#include "utils.hpp"


static bool isInt(char* c)
{
    for (int i=0; c[i]; i++)
    {
        if (!isdigit(c[i]))
            return false;
    }
    return true;
}


static bool parsing(int argc, char* argv[], Config& cf)
{
    const struct option opts[] = {
        {"i", 1, NULL, 'i'},
        {"o", 1, NULL, 'o'},
        {"p", 1, NULL, 'p'},
        {"g", 1, NULL, 'g'},
        {"d", 1, NULL, 'd'},
        {"h", 1, NULL, 'h'},
        {"r", 1, NULL, 'r'},
        {"opt", 1, NULL, 'O'}
    };
    const char *optstring = "i:o:pgdhO";
    int c;

    while ((c = getopt_long(argc, argv, optstring, opts, NULL)) != -1)
    {
        switch (c)
        {
            case 'i':
                // if (!access(optarg, R_OK))
                // {
                //     printf ("-i: Input file is unreadable or does not exist\n");
                //     return false;
                // }
                cf.infile = optarg;
                break;
            case 'o':
                // ofstream out(optarg);
                // if (out.fail())
                // {
                //     printf ("-o: Output file opening failed\n");
                //     return false;
                // }
                // out.close();
                cf.outfile = optarg;
                break;
            case 'p':
                if (!isInt(optarg))
                {
                    printf ("--p: Please enter a positive integer\n");
                    return false;
                }
                cf.Pc = atoi(optarg);
                break;
            case 'g':
                if (!isInt(optarg))
                {
                    printf ("--g: Please enter a positive integer\n");
                    return false;
                }
                cf.generation = atoi(optarg);
                break;
            case 'd':
                if (!isInt(optarg))
                {
                    printf ("--d: Please enter a positive integer\n");
                    return false;
                }
                cf.max_descendant = atoi(optarg);
                break;
            case 'h':
                if (!isInt(optarg))
                {
                    printf ("--h: Please enter a positive integer\n");
                    return false;
                }
                cf.ymin = atoi(optarg);
                break;
            case 'r':
                if (!isInt(optarg))
                {
                    printf ("--r: Please enter a positive integer\n");
                    return false;
                }
                cf.rank = atoi(optarg);
                break;
            case 'O':
                if (!isInt(optarg))
                {
                    printf ("--opt: Please enter a positive integer\n");
                    return false;
                }
                cf.to_optimize = atoi(optarg);
                break;
            case '?':
                printf("Unknown option %s\n", optarg);
                break;
            default:
                printf("------------------------------------\n");
        }
    }

    return true;
}