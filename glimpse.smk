



# GLIMPSE has no working conda or anything
# to make installation backwards compatible, download static binaries from release
rule install_glimpse:
    output:
        "glimpse/GLIMPSE_chunk_static",
        "glimpse/GLIMPSE_phase_static",
    shell:
        "wget -O glimpse/GLIMPSE_chunk_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static && "
        "wget -O glimpse/GLIMPSE_phase_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static && "
        "chmod a+x glimpse/GLIMPSE_*_static "


