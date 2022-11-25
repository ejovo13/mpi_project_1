#include "cli.h"

args_t *new_args(int argc, char **argv) {

    args_t *arg = (args_t *) malloc(sizeof(*arg));

    arg->argc = argc;
    arg->argv = argv;
    arg->lmodel = -1;
    arg->lbin = -1;
    arg->size_dataset = NULL;
    arg->data = NULL;
    arg->from = false;
    arg->a = 0;
    arg->txt = false;
    arg->predict = false;
    arg->diff = false;
    arg->recompute = false;

    return arg;
}

void usage(char ** argv)
{
    printf("%s [OPTIONS]\n\n", argv[0]);
    printf("Options:\n");
    printf("--[small | med | hi | ultra]    dataset to model\n");
    printf("--lmodel l                      degree of the model to compute\n");
    printf("--lbin L                        degree of the stored binary file\n");
    printf("[--from a]                      compute a model of degree lmodel starting from degree a\n");
    printf("                                the difference between the model and predicted value\n");
    printf("[--txt]                         output the model as a text file\n");
    printf("[--predict]                     predict the altitude values (f_hat)\n");
    printf("[--diff]                        predict the altitude values (f_hat) and compute\n");
    printf("[--recompute]                   compute the model coefficients no matter what files\n");
    printf("                                are present\n");
    printf("\n");
    exit(0);
}

void print_args(const args_t *args) {
    printf("args = {\n");
    printf("  .size_dataset\t = %s\n", args->size_dataset);
    printf("  .lmodel\t = %d\n", args->lmodel);
    printf("  .lbin\t\t = %d\n", args->lmodel);
    printf("  .from\t\t = %s\n", btos(args->from));
    printf("  .a\t\t = %d\n", args->a);
    printf("  .txt\t\t = %s\n", btos(args->txt));
    printf("  .predict\t = %s\n", btos(args->predict));
    printf("  .diff\t\t = %s\n", btos(args->diff));
    printf("  .recompute\t = %s\n", btos(args->recompute));
    printf("}\n");

}

args_t *process_command_line_options(int argc, char ** argv, bool __log)
{
    struct option longopts[] = {
        {"lmodel", required_argument, NULL, 'l'},
        {"lbin", required_argument, NULL, 'L'},
        {"small", no_argument, NULL, 's'},
        {"med", no_argument, NULL, 'm'},
        {"hi", no_argument, NULL, 'h'},
        {"ultra", no_argument, NULL, 'u'},
        {"txt", no_argument, NULL, 't'},
        {"predict", no_argument, NULL, 'p'},
        {"diff", no_argument, NULL, 'd'},
        {"from", required_argument, NULL, 'f'},
        {"recompute", no_argument, NULL, 'r'},
        {NULL, 0, NULL, 0}
    };

    args_t *args = new_args(argc, argv);

    char ch;
    while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (ch) {

        case 'L':
            args->lbin = atoll(optarg);
            break;
        case 'l':
            args->lmodel = atoll(optarg);
            args->lbin = args->lmodel;
            break;
        case 's':
            args->size_dataset = "small";
            args->data = get_data_small(__log);
            break;
        case 'm':
            args->size_dataset = "med";
            args->data = get_data_med(__log);
            break;
        case 'h':
            args->size_dataset = "hi";
            args->data = get_data_hi(__log);
            break;
        case 'u':
            args->size_dataset = "ultra";
            args->data = get_data_ultra(__log);
            break;
        case 't':
            args->txt = true;
            break;
        case 'p':
            args->predict = true;
            break;
        case 'd':
            args->diff = true;
            args->predict = true;
            break;
        case 'f':
            args->from = true;
            args->a = atoi(optarg);
            break;
        case 'r':
            args->recompute = true;
        break;
        default:
            errx(1, "Unknown option\n");
        }
    }

    /* missing required args? */
    if (args->size_dataset == NULL || args->data == NULL || args->lbin < 0 || args->lmodel < 0)
        usage(argv);

    return args;
}