#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <signal.h>
#include <string.h>
#include <errno.h>

int diag_size = 0;

void diag_reset(void) {
    /* TODO: find a portable way to do this without looping */
    while (diag_size) {
        printf("\x08 \x08");
        --diag_size;
    }
}

void keep_diag(void) {
    if (diag_size)
        (void) !write(1, "\n", 1);
    diag_size = 0;
}

unsigned int diag(char *format, ...) {
    diag_reset();

    va_list ap;
    va_start(ap, format);
    diag_size += vprintf(format, ap);
    va_end(ap);
    fflush(stdout);
    return diag_size;
}

int diag_fail(char *msg) {
    (void) !write(2, msg, strlen(msg));
    return 0;
}

/* needs to be async-signal-safe */
void diag_fatal(char *msg) {
    (void) !write(2, msg, strlen(msg));
    exit(1);
}

/* Adapted from:
 *   https://man7.org/tlpi/code/online/dist/pgsjc/handling_SIGTSTP.c.html
 * This code is copyright 2022, Michael Kerrisk, and is licensed under the
 * GNU General Public License, version 3.
 */
void diag_TSTP(int sig) {
    sigset_t tstpMask, prevMask;
    int savedErrno;
    struct sigaction sa;

    savedErrno = errno;                 /* In case we change 'errno' here */

    /* no need to keep_diag(), we'll get a newline for free */
    diag_size = 0;

    if (signal(SIGTSTP, SIG_DFL) == SIG_ERR)
        diag_fatal("signal");           /* Set handling to default */

    raise(SIGTSTP);                     /* Generate a further SIGTSTP */

    /* Unblock SIGTSTP; the pending SIGTSTP immediately suspends the program */
    sigemptyset(&tstpMask);
    sigaddset(&tstpMask, SIGTSTP);
    if (sigprocmask(SIG_UNBLOCK, &tstpMask, &prevMask) == -1)
        diag_fatal("sigprocmask");

    /* Execution resumes here after SIGCONT */
    if (sigprocmask(SIG_SETMASK, &prevMask, NULL) == -1)
        diag_fatal("sigprocmask");      /* Reblock SIGTSTP */

    sigemptyset(&sa.sa_mask);           /* Reestablish handler */
    sa.sa_flags = SA_RESTART;
    sa.sa_handler = diag_TSTP;
    if (sigaction(SIGTSTP, &sa, NULL) == -1)
        diag_fatal("sigaction");

    errno = savedErrno;
}

int init_diag(void) {
    struct sigaction sa;

    /* Only establish handler for SIGTSTP if it is not being ignored */
    if (sigaction(SIGTSTP, NULL, &sa) == -1)
        return diag_fail("sigaction fetch");

    if (sa.sa_handler != SIG_IGN) {
        sigemptyset(&sa.sa_mask);
        sa.sa_flags = SA_RESTART;
        sa.sa_handler = diag_TSTP;
        if (sigaction(SIGTSTP, &sa, NULL) == -1)
            return diag_fail("sigaction put");
    }
    return 1;
}

