
use htslib::bcf;

use call::site::Site;
use Prob;

pub fn write(writer: &mut bcf::Writer, mut calls: Vec<(Site, Prob)>) {
    for (site, prob) in calls.drain() {
        let record = site.into_record(prob);
        writer.write(&record).ok().expect("Error writing calls.");
    }
}
