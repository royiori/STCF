  Universe             = vanilla
  Notification         = Never
  GetEnv               = True
        Executable           = EXECPATH
        Arguments            = ARGU1 ARGU2
  Output               = OUT
  Error                = ERR
  Log                  = LOG
  +Group               = "BESIII"
  should_transfer_files= yes
        requirements         = (substr(Machine,0,4)=="bl-3"||substr(Machine,0,4)=="bl-4"||substr(Machine,0,4)=="bl-5")&&(machine != "bl-2-12.hep.ustc.edu.cn")
  WhenToTransferOutput = ON_EXIT_OR_EVICT
        accounting_group                = long
  OnExitRemove         = TRUE
  Queue
