#include "ittnotify.h"

void fortran_itt_resume()
{
    __itt_resume();
}

void fortran_itt_pause()
{
    __itt_pause();
}

void fortran_itt_detach()
{
    __itt_detach();
}

__itt_domain* fortran_itt_domain_create(const char *name)
{
    __itt_domain *newdomain = __itt_domain_create(name);
    newdomain->flags = 1; /* enable the domain */
    return newdomain;
}

__itt_string_handle* fortran_itt_string_handle_create(const char *name)
{
    __itt_string_handle* ptr_string = __itt_string_handle_create(name);
    return ptr_string;
}

/* I didn't understand what the intended use of the ID is so far, thus I
  skipped this extra argument, the official examples also use NULL.*/
void fortran_itt_frame_begin(const __itt_domain *domain)
{
    __itt_frame_begin_v3(domain, NULL);
}

void fortran_itt_frame_end(const __itt_domain *domain)
{
    __itt_frame_end_v3(domain, NULL);
}

void fortran_itt_task_begin(const __itt_domain *domain, __itt_string_handle *name)
{
    __itt_task_begin(domain, __itt_null, __itt_null, name);
}

void fortran_itt_task_end(const __itt_domain *domain)
{
    __itt_task_end(domain);
}
