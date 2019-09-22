## ---- prp-setup ---------------------
prp_base <- function(path.prp, prefix) {
    df.prp.sw <- ldply(c("cdnc", "cf", "lwp"), function(pert) {
        fname <- sprintf("%s/%s_flsh_diff_%s_timmean.nc", path.prp, prefix, pert)
        varname <- sprintf("flsh_diff_%s", pert)
        nc <- nc_open(fname)
        on.exit(nc_close(nc))
        var <- ncvar_get(nc, varname)
        lon <- ncvar_get(nc, "lon")
        lat <- ncvar_get(nc, "lat")
        expand.grid(lon = as.vector(lon), lat = as.vector(lat)) %>%
            mutate(erf = -as.vector(var),
                   pert = pert,
                   spectrum = "SW")
    })

    df.prp.lw <- ldply(c("cdnc", "cf", "lwp"), function(pert) {
        fname <- sprintf("%s/%s_flth_diff_%s_timmean.nc", path.prp, prefix, pert)
        varname <- sprintf("flth_diff_%s", pert)
        nc <- nc_open(fname)
        on.exit(nc_close(nc))
        var <- ncvar_get(nc, varname)
        lon <- ncvar_get(nc, "lon")
        lat <- ncvar_get(nc, "lat")
        expand.grid(lon = as.vector(lon), lat = as.vector(lat)) %>%
            mutate(erf = -as.vector(var),
                   pert = pert,
                   spectrum = "LW")
    })

    bind_rows(df.prp.sw, df.prp.lw)
}

prp <- function(path.prp) {
    prp_base(path.prp, "prp")
}

prp.fw <- function(path.prp) {
    prp_base(path.prp, "prp_fw")
}

prp.bk <- function(path.prp) {
    prp_base(path.prp, "prp_bk")
}

prp.global <- function(path.prp) {
    df.prp.sw <- ldply(c("cdnc", "cf", "lwp"), function(pert) {
        fname <- sprintf("%s/prp_flsh_diff_%s_global.nc", path.prp, pert)
        varname <- sprintf("flsh_diff_%s", pert)
        nc <- nc_open(fname)
        on.exit(nc_close(nc))
        var <- ncvar_get(nc, varname)
        data.frame(erf = -as.vector(var)) %>%
            mutate(pert = pert,
                   spectrum = "SW")
    })
    
    df.prp.lw <- ldply(c("cdnc", "cf", "lwp"), function(pert) {
        fname <- sprintf("%s/prp_flth_diff_%s_global.nc", path.prp, pert)
        varname <- sprintf("flth_diff_%s", pert)
        nc <- nc_open(fname)
        on.exit(nc_close(nc))
        var <- ncvar_get(nc, varname)
        data.frame(erf = -as.vector(var)) %>%
            mutate(pert = pert,
                   spectrum = "LW")
    })
    
    bind_rows(df.prp.sw, df.prp.lw)
}

prp.combine.lw.and.sw <- function(df) {
    df %>%
        tidyr::spread(spectrum, erf) %>%
        dplyr::mutate(erf = LW + SW) %>%
        dplyr::select(-c(LW, SW))
}

prp.plot <- function(df, range = 5, palette = "RdYlBu", symmetric = TRUE, direction = -1, title = "W~m$^{-2}$") {
    if (length(range) == 1) {
        range <- c(-range,
                   ifelse(symmetric, range, 0))
    }
    df %>%
        mutate(lon = ifelse(lon <= 180, lon, lon - 360)) %>%
        ggplot(aes(lon, lat, fill = pmax(pmin(erf, range[2]), range[1]))) +
        geom_raster() +
        scale_x_geo(facet = FALSE) + scale_y_geo() +
        coord_fixed(xlim = c(-180, 180), ylim = c(-80, 80), expand = FALSE) +
        ## scale_x_continuous("", labels = NULL, breaks = NULL) +
        ## scale_y_continuous("", labels = NULL, breaks = NULL) +
        ## scale_fill_manual(values = col.frac, name = expression(f[liq])) +
        ## scale_fill_warmfrac() +
        ## scale_fill_brewer("$F_\\mathcal{L}~(\\text{W~m}^{-2})$", palette = "RdBu", drop = FALSE, direction = -1) +
        scale_fill_distiller(labels = tikz_sanitize, title, palette = palette, 
                             limits = range, direction = direction) +
        geom_world_polygon(highres = FALSE) +
        theme_bw(12) +
        facet_grid(as.character(dplyr::groups(df))[1:2] %>%
                   coalesce(".") %>%
                   Reduce(function(x, y) paste(x, y, sep = " ~ "), .)) +
        theme(legend.position = "bottom", legend.box = "horizontal") +
        guides(fill = guide_colorbar(direction = "horizontal", title.vjust = 0.75, barwidth = 8))
}

prp.contour <- function(df, range = 5, binwidth = 0.5, palette = "RdYlBu", symmetric = TRUE, direction = -1, title = "W~m$^{-2}$") {
    if (length(range) == 1) {
        range <- c(-range,
                   ifelse(symmetric, range, 0))
    }
    df %>%
        mutate(lon = ifelse(lon <= 180, lon, lon - 360)) %>%
        ggplot(aes(lon, lat)) +
        geom_contour_fill(aes(z = pmax(pmin(erf, range[2]), range[1])),
                          binwidth = binwidth, na.fill = TRUE) +
        scale_x_geo(facet = FALSE) + scale_y_geo() +
        coord_fixed(xlim = c(-180, 180), ylim = c(-80, 80), expand = FALSE) +
        ## scale_x_continuous("", labels = NULL, breaks = NULL) +
        ## scale_y_continuous("", labels = NULL, breaks = NULL) +
        ## scale_fill_manual(values = col.frac, name = expression(f[liq])) +
        ## scale_fill_warmfrac() +
        ## scale_fill_brewer("$F_\\mathcal{L}~(\\text{W~m}^{-2})$", palette = "RdBu", drop = FALSE, direction = -1) +
        scale_fill_distiller(labels = tikz_sanitize, title, palette = palette,
                             limits = range, direction = direction,
                             breaks = MakeBreaks(binwidth),
                             guide = guide_colorstrip(title.vjust = 0.75, title.hjust = 0, barwidth = 8)) +
        geom_world_polygon(highres = FALSE) +
        theme_bw(12) +
        facet_grid(as.character(dplyr::groups(df))[1:2] %>%
                   coalesce(".") %>%
                   Reduce(function(x, y) paste(x, y, sep = " ~ "), .)) +
        theme(legend.position = "bottom", legend.box = "horizontal") 
}

delta.ln.nd <- function(path.prp) {
    df.delta.ln.nd <-
        ldply(c("CDNC_BURDEN"), function(varname) {
            fname.pd <- sprintf("%s/cdnc_burden_pd_timmean.nc", path.prp)
            fname.pi <- sprintf("%s/cdnc_burden_pi_timmean.nc", path.prp)
            nc.pd <- nc_open(fname.pd)
            nc.pi <- nc_open(fname.pi)
            on.exit(nc_close(nc.pd))
            on.exit(nc_close(nc.pi))
            pd <- ncvar_get(nc.pd, varname)
            pi <- ncvar_get(nc.pi, varname)
            lon <- ncvar_get(nc.pd, "lon")
            lat <- ncvar_get(nc.pd, "lat")
            expand.grid(lon = as.vector(lon), lat = as.vector(lat)) %>%
                mutate(delta.ln.nd = as.vector(log(pd/pi)))
    })
}

clouds.pd <- function(path.prp) {
    fname.pd <- sprintf("%s/cdnc_burden_pd_timmean.nc", path.prp)
    nc.pd <- nc_open(fname.pd)
    on.exit(nc_close(nc.pd))
    nd <- ncvar_get(nc.pd, "CDNC_BURDEN")
    lwp <- ncvar_get(nc.pd, "xlvi")
    fc <- ncvar_get(nc.pd, "aclcov")
    lon <- ncvar_get(nc.pd, "lon")
    lat <- ncvar_get(nc.pd, "lat")
        expand.grid(lon = as.vector(lon), lat = as.vector(lat)) %>%
            mutate(nd = as.vector(nd),
                   lwp = as.vector(lwp),
                   fc = as.vector(fc))
}

clouds.pi <- function(path.prp) {
    fname.pi <- sprintf("%s/cdnc_burden_pi_timmean.nc", path.prp)
    nc.pi <- nc_open(fname.pi)
    on.exit(nc_close(nc.pi))
    nd <- ncvar_get(nc.pi, "CDNC_BURDEN")
    lwp <- ncvar_get(nc.pi, "xlvi")
    fc <- ncvar_get(nc.pi, "aclcov")
    lon <- ncvar_get(nc.pi, "lon")
    lat <- ncvar_get(nc.pi, "lat")
        expand.grid(lon = as.vector(lon), lat = as.vector(lat)) %>%
            mutate(nd = as.vector(nd),
                   lwp = as.vector(lwp),
                   fc = as.vector(fc))
}

kernel.sensitities <- function(fname) {
    nc <- try(ncdf4::nc_open(fname), silent = TRUE)
    if (class(nc) == "try-error")
        return(NULL)
    t <- ncdf4::ncvar_get(nc, "time")
    lon <- ncdf4::ncvar_get(nc, "lon")
    lat <- ncdf4::ncvar_get(nc, "lat")
    xlvi <- ncdf4::ncvar_get(nc, "XLVI")
    fsw_diff <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_LWP")
    fsw_total_top_unpert <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_UNPERT")
    flw_total_top_unpert <- ncdf4::ncvar_get(nc, "FLW_TOTAL_TOP_UNPERT")
    fsw_total_top_lwp    <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_LWP")
    flw_total_top_lwp    <- ncdf4::ncvar_get(nc, "FLW_TOTAL_TOP_LWP")
    fsw_total_top_cdnc   <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_CDNC")
    flw_total_top_cdnc   <- ncdf4::ncvar_get(nc, "FLW_TOTAL_TOP_CDNC")
    fsw_total_top_cf     <- ncdf4::ncvar_get(nc, "FSW_TOTAL_TOP_CLDFRA")
    flw_total_top_cf     <- ncdf4::ncvar_get(nc, "FLW_TOTAL_TOP_CLDFRA")
    cos_mu0 <- ncdf4::ncvar_get(nc, "COS_MU0")
    xlvi       <- ncdf4::ncvar_get(nc, "XLVI")
    cdnc       <- ncdf4::ncvar_get(nc, "CDNC")
    cldfra     <- ncdf4::ncvar_get(nc, "CLDFRA")
    cldfra_liq <- ncdf4::ncvar_get(nc, "CLDFRA_LIQ")
    ncdf4::nc_close(nc)
    ## finite difference approximation to logarithmic derivatives
    dfsw.dlog.lwp <-
        2 * (fsw_total_top_lwp - fsw_total_top_unpert) /
        (fsw_total_top_lwp + fsw_total_top_unpert) /
        (0.1 / 1.05)
    dfsw.dlog.cdnc <-
        2 * (fsw_total_top_cdnc - fsw_total_top_unpert) /
        (fsw_total_top_cdnc + fsw_total_top_unpert) /
        (0.1 / 1.05)
    dfsw.dcf <-
        2 * (fsw_total_top_cdnc - fsw_total_top_unpert) /
        (fsw_total_top_cdnc + fsw_total_top_unpert) /
        (0.1 * cldfra_liq)
    expand.grid(lon = as.vector(lon),
                lat = as.vector(lat),
                time = as.vector(t)) %>%
        dplyr::mutate(dfsw.dlog.lwp  = as.vector(dfsw.dlog.lwp ),
                      dfsw.dlog.cdnc = as.vector(dfsw.dlog.cdnc),
                      dfsw.dcf       = as.vector(dfsw.dcf      ),
                      fsw_total_top_unpert = as.vector(fsw_total_top_unpert),
                      flw_total_top_unpert = as.vector(flw_total_top_unpert),
                      cos_mu0            = as.vector(cos_mu0           ),
                      xlvi               = as.vector(xlvi              ),
                      cdnc               = as.vector(cdnc              ),
                      cldfra             = as.vector(cldfra            ),
                      cldfra_liq         = as.vector(cldfra_liq        ))
}

gamma <- 2^(1/2) * 3^(4/3) / 5 * pi^(-1/3) * (1e3)^(-2/3) / 1.1 * (1e-6)^-(1/6)

tau <- function(N, L) {
    gamma * N^(1/3) * L^(5/6)
}

alpha <- function(tau, a1 = 0.092, a2 = 1.43, g = 0.85) {
    (tau * (1 - g) + a1) / (tau * (1 - g) + a2)
}

alpha.map <- function(N, L) {
    expand.grid(N = N, L = L) %>%
        ## calculate tau and alpha on the N, L grid
        mutate(tau = tau(N, L),
               alpha = alpha(tau)) %>%
        ## calculate dalpha/dN and dalpha/dL on the N, L grid
        group_by(N) %>%
        mutate(dalpha.dL = (alpha - lag(alpha)) / (L - lag(L)),
               dalpha.dlogL = (alpha - lag(alpha)) / (log(L) - log(lag(L)))) %>%
        ungroup() %>%
        group_by(L) %>%
        mutate(dalpha.dN = (alpha - lag(alpha)) / (N - lag(N)),
               dalpha.dlogN = (alpha - lag(alpha)) / (log(N) - log(lag(N)))) %>%
        ungroup() %>%
        filter(!is.na(dalpha.dL),
               !is.na(dalpha.dN))

}
