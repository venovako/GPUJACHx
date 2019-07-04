#ifndef USE_STRAT_ARRAY_DECLARATOR
#define USE_STRAT_ARRAY_DECLARATOR
#endif /* !USE_STRAT_ARRAY_DECLARATOR */
#ifndef EXPORT_VAR
#define EXPORT_VAR __attribute__((visibility("default")))
#else /* EXPORT_VAR */
#error EXPORT_VAR already defined
#endif /* ?EXPORT_VAR */
