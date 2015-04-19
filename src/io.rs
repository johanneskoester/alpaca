
use htslib::bcf;

use call::site::Site;
use LogProb;

pub fn write(writer: &mut bcf::Writer, mut calls: Vec<(Site, LogProb)>) {
    for (site, prob) in calls.drain() {
        let record = site.into_record(prob, &writer.header);
        writer.write(&record).ok().expect("Error writing calls.");
    }
}
