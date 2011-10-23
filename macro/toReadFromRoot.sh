#!/bin/zsh

writeScript() {

    root -b <<EOF

.L makeMCPlots.C++
makeMCPlots("Ele",1545,true,1)
.q
EOF
    
}

writeScript;




