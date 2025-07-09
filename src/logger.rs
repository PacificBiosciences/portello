use crate::globals::PROGRAM_NAME;

/// If debug is true set the default logger to the more verbose debug level
///
pub fn setup_logger(debug: bool) -> Result<(), fern::InitError> {
    let level = if debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    let logger = fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{}[{}][{}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                PROGRAM_NAME,
                record.level(),
                message
            ))
        })
        .level(level)
        .chain(std::io::stderr());

    logger.apply()?;
    Ok(())
}
